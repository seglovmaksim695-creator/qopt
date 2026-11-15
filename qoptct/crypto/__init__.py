# SPDX-License-Identifier: Apache-2.0
"""qoptct.crypto: building blocks of cryptographic circuits on logical ANDs (Gidney 1709.06648: AND = 4 T, erased by
X-measurement with a Clifford correction, 0 T), each checked exactly at build time.

Modules:
    lookup   narrow table lookups z ^= T[a] (bilinear split, QROAM s=1, ANF/Davio), measurement-based unlookup,
             lookup-accumulate acc += T[a]
"""
