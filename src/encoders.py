"""Encoder registry. Loading code, preprocessing, and output convention follow each HF model card.

Each loader returns (model, transform, forward, dim, info):
  transform: PIL.Image (RGB) -> float tensor [3, 224, 224]
  forward:   model, batch tensor -> [B, dim] embedding
"""
import os

import torch

# Pinned HF revisions (commit sha) used for every extraction; recorded in emb/<encoder>/meta.json.
REVISIONS = {
    "owkin/phikon": "057cc0295895c2df3dd7681a89680da6015cbefe",
    "owkin/phikon-v2": "2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd",
    "MahmoodLab/UNI": "b55a5ec6cade1a39edfe6534189a9b8ca7a022f0",
    "MahmoodLab/UNI2-h": "d517a8dd47902dd7c308b3c36f63bce47e7b9a43",
    "paige-ai/Virchow2": "3158645804b69e3f3bc4439d4116edddf0840a72",
    "prov-gigapath/prov-gigapath": "64f9e26c15019f2d4f6d9113c6822f88bb16b01b",
}


def _hub(repo):
    """timm hub id pinned to the revision (timm splits 'hf-hub:org/name@rev')."""
    return f"hf-hub:{repo}@{REVISIONS[repo]}"


class HFProcessorTransform:
    """Picklable wrapper (DataLoader workers are spawned on macOS)."""

    def __init__(self, proc):
        self.proc = proc

    def __call__(self, img):
        return self.proc(img, return_tensors="pt")["pixel_values"][0]


def _hf_cls_forward(m, x):
    return m(pixel_values=x).last_hidden_state[:, 0, :]  # CLS token, per model card


def _hf_transformers(repo):
    from transformers import AutoImageProcessor, AutoModel

    proc = AutoImageProcessor.from_pretrained(repo, revision=REVISIONS[repo])
    model = AutoModel.from_pretrained(repo, revision=REVISIONS[repo]).eval()
    return model, HFProcessorTransform(proc), _hf_cls_forward, model.config.hidden_size, {"processor": proc.to_dict()}


def load_phikon():
    out = _hf_transformers("owkin/phikon")
    return (*out[:4], {**out[4], "repo": "owkin/phikon", "output": "last_hidden_state[:,0] (CLS)"})


def load_phikon_v2():
    out = _hf_transformers("owkin/phikon-v2")
    return (*out[:4], {**out[4], "repo": "owkin/phikon-v2", "output": "last_hidden_state[:,0] (CLS)"})


def _timm_cfg_transform(model):
    from timm.data import resolve_data_config
    from timm.data.transforms_factory import create_transform

    cfg = resolve_data_config(model.pretrained_cfg, model=model)
    return create_transform(**cfg), {k: (list(v) if isinstance(v, tuple) else v) for k, v in cfg.items()}


def _cls_forward(m, x):
    return m(x)  # timm num_classes=0, global_pool='token' -> CLS


def load_uni():
    import timm

    model = timm.create_model(_hub("MahmoodLab/UNI"), pretrained=True, init_values=1e-5, dynamic_img_size=True).eval()
    tf, cfg = _timm_cfg_transform(model)
    return model, tf, _cls_forward, 1024, {"repo": "MahmoodLab/UNI", "timm_data_cfg": cfg, "transform": str(tf),
                                          "output": "timm forward, global_pool=token (CLS)"}


def load_uni2_h():
    import timm

    kw = dict(img_size=224, patch_size=14, depth=24, num_heads=24, init_values=1e-5, embed_dim=1536,
              mlp_ratio=2.66667 * 2, num_classes=0, no_embed_class=True, mlp_layer=timm.layers.SwiGLUPacked,
              act_layer=torch.nn.SiLU, reg_tokens=8, dynamic_img_size=True)
    model = timm.create_model(_hub("MahmoodLab/UNI2-h"), pretrained=True, **kw).eval()
    tf, cfg = _timm_cfg_transform(model)
    return model, tf, _cls_forward, 1536, {"repo": "MahmoodLab/UNI2-h", "timm_data_cfg": cfg, "transform": str(tf),
                                          "output": "timm forward, global_pool=token (CLS)"}


def load_virchow2():
    import timm
    from timm.layers import SwiGLUPacked

    model = timm.create_model(_hub("paige-ai/Virchow2"), pretrained=True, mlp_layer=SwiGLUPacked,
                              act_layer=torch.nn.SiLU).eval()
    tf, cfg = _timm_cfg_transform(model)

    def forward(m, x):
        out = m(x)  # [B, 261, 1280]: CLS, 4 register tokens, 256 patch tokens
        return torch.cat([out[:, 0], out[:, 5:].mean(1)], dim=-1)

    return model, tf, forward, 2560, {"repo": "paige-ai/Virchow2", "timm_data_cfg": cfg, "transform": str(tf),
                                     "output": "concat(CLS, mean(patch tokens[5:])), 2560-d"}


def load_gigapath():
    import timm
    from torchvision import transforms

    model = timm.create_model(_hub("prov-gigapath/prov-gigapath"), pretrained=True).eval()
    tf = transforms.Compose([
        transforms.Resize(256, interpolation=transforms.InterpolationMode.BICUBIC),
        transforms.CenterCrop(224),
        transforms.ToTensor(),
        transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)),
    ])
    return model, tf, _cls_forward, 1536, {"repo": "prov-gigapath/prov-gigapath", "transform": str(tf),
                                          "output": "tile encoder timm forward (CLS), 1536-d"}


FLASH_REPO, FLASH_REV = "prov-gigapath/prov-gigapath-flash", "f5c3a709ad81d3a5605bb2eb38cd5c2e53da0375"
FLASH_SHA256 = "ea5e36c2b3593b073868f91b9016559a620f0d6752e3861cb2bb0b5312238473"
FLASH_BATCH = 128  # every forward runs at this fixed batch shape (padded), so fp16 outputs do not depend on batching


def load_flash():
    """Prov-GigaPath-flash tile encoder (public, gated HF weights; ViT-S/16 DINOv2-small with SwiGLU and LayerScale,
    21.7 M params). Numerics: uint8 in, fp32 ImageNet normalisation, fp16 weights and input, fixed batch shape;
    output CLS || mean(patch tokens) = 768-d, pooled in fp32."""
    import hashlib

    from huggingface_hub import hf_hub_download
    from timm.layers import SwiGLUPacked
    from timm.models.vision_transformer import VisionTransformer
    from torchvision.transforms import PILToTensor

    path = hf_hub_download(FLASH_REPO, "pytorch_model.bin", revision=FLASH_REV)
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()
    assert h == FLASH_SHA256, f"unexpected prov-gigapath-flash weights {h[:12]}"
    model = VisionTransformer(img_size=224, patch_size=16, embed_dim=384, depth=12, num_heads=6,
                              mlp_ratio=2048 / 384.0, init_values=1e-5, num_classes=0, global_pool="token",
                              class_token=True, reg_tokens=0, mlp_layer=SwiGLUPacked, act_layer=torch.nn.SiLU)
    sd = torch.load(path, map_location="cpu", weights_only=True)
    model.load_state_dict(sd.get("model", sd) if "patch_embed.proj.weight" not in sd else sd, strict=True)
    model = model.eval().half()
    mean = torch.tensor((0.485, 0.456, 0.406)).view(1, 3, 1, 1)
    std = torch.tensor((0.229, 0.224, 0.225)).view(1, 3, 1, 1)

    def forward(m, x):  # x: uint8 [B,3,224,224] on device
        out = []
        mu, sd_ = mean.to(x.device), std.to(x.device)
        for s in range(0, len(x), FLASH_BATCH):
            chunk = x[s:s + FLASH_BATCH]
            k = len(chunk)
            if k < FLASH_BATCH:
                chunk = torch.cat([chunk, chunk.new_zeros((FLASH_BATCH - k, *chunk.shape[1:]))])
            z = ((chunk.float().div_(255.0) - mu) / sd_).half()
            f = m.forward_features(z)[:k]
            out.append(torch.cat([f[:, 0].float(), f[:, 1:].float().mean(1)], -1))
        return torch.cat(out)

    return model, PILToTensor(), forward, 768, {
        "repo": FLASH_REPO, "weights_sha256": FLASH_SHA256,
        "output": "CLS || mean(patch tokens), 768-d",
        "numerics": f"fp16 weights, fp32 normalise, fp16 input, fixed batch {FLASH_BATCH}"}


ENCODERS = {
    "phikon": load_phikon,
    "phikon_v2": load_phikon_v2,
    "uni": load_uni,
    "uni2_h": load_uni2_h,
    "virchow2": load_virchow2,
    "gigapath": load_gigapath,
    "flash": load_flash,
}


def load(name):
    model, tf, fwd, dim, info = ENCODERS[name]()
    info["revision"] = REVISIONS.get(info["repo"], FLASH_REV if info["repo"] == FLASH_REPO else None)
    return model, tf, fwd, dim, info
