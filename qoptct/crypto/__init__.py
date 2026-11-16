# Copyright 2026 Maksym Shchehlov. Licensed under the Apache License, Version 2.0.
"""qoptct.crypto -- verified arithmetic circuit generators for cryptographic fields (qoptct v0.5).

Currently: GF(2^m) (binary fields of the NIST/SEC curves), see qoptct.crypto.gf2m."""
from .gf2m import (default_modulus, gf2m_bilinear, gf2m_multiply, gf2m_mult_count, gf2m_square, gf2m_inverse,
                   gf2m_divide, gf2m_inverse_count, gf2m_divide_count, gf2m_itoh_gate_check, verify_multiplier)

__all__ = ["default_modulus", "gf2m_bilinear", "gf2m_multiply", "gf2m_mult_count", "gf2m_square", "gf2m_inverse",
           "gf2m_divide", "gf2m_inverse_count", "gf2m_divide_count", "gf2m_itoh_gate_check", "verify_multiplier"]
