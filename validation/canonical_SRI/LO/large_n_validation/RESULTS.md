# LO high-order validation

All checks below passed in this execution.

| Check | Range | Arithmetic |
| --- | --- | --- |
| Archived file integrity | 11 original files | SHA-256 |
| SRI / current Bell / archived Bell / Riccati | n=0,...,8 (9 expressions) | symbolic, zero residual |
| Current / archived Bell / Riccati and direct partial Bell | n=0,...,200; 201 values; direct sum through n=25 or requested maximum | exact rational |
| Current / archived Bell / Riccati | n=0,...,1000; 9 cases across three primes | exact finite-field |
| Recovered projective jet / Riccati | n=0,...,1000; T_1,...,T_1001 | exact rational; seed 4321 |

Elapsed time: 54.74 seconds.

`T_{n+1}=alpha_n^t`, `S_j=alpha_{j-1}^s`, and `T_1=-r/8`.

The n=1000 stress tests concern the LO derivative hierarchy, not SR order p=1000.
Finite-field comparisons are exact in each field at the stated inputs; they do not
replace a symbolic polynomial identity. SRI itself is compared symbolically through n=8.

## Source fingerprints

```text
0ec97d8fe443c934868716872a4b250c990f96178a6553c9d77bcf96913de2c2  source/__init__.py
2ef9747e933aaaed6ffb04c9b8d6494c575836792d5a4d05fb79279cea37623a  source/__main__.py
ca6844dd1b97a9a43fe604c8c0df7695dd84ebfe9ba6d9e41117bc6601d3bc11  source/background.py
08974a75836fb08d52a2b446409fa3aba3c262db6415148bdb0fde68a93aaada  source/cli.py
5501e45bb679064f91834fde800c0ab3f338657876b56e5badc8706e61eefbd0  source/closed/__init__.py
7b6e9325a47c465bd5611e4b81b90fd810573105ea8dda0e28a43b5d653f06c5  source/closed/algebra.py
4ee5c3fb42fbac637c89d4eae63ec35fd64c402a8875c158cc587e08b86c506f  source/closed/consistency.py
3dd4020c88f39d9f0969bd5ac11bf7f271738d9ab7936ad33f0a25b367b56f1b  source/closed/inversion.py
77c84e05e0d18e19edd1a8e3a98c3c8121bb5104de6d30df373e9529bc6abfa8  source/closed/spectra.py
8c61d56d826831fbb4b5a148edaceb32ba593b057795afd2aecd3532b789ccb7  source/consistency.py
0e2cfd8712f57c4186a920071733254d7fae0285c521365999acfc45644a154c  source/generic_modes.py
9ded71d9731126cbaafef6d7b952be677a5a357982cd625714cd8bc65ab2749a  source/inversion.py
07c7c626ecdfe0b1296d244f03c94890b16be592c1f23878e74c99897c5730bb  source/iterative/__init__.py
efc9f0d571ffaddc0032b88e034ed10731df31cc48f1f9d19e296d4fa5e1ea8d  source/iterative/background.py
2e924d76666ca8665869efae73f4947de2663dbdaacb916bc486c35647d9d4eb  source/iterative/consistency.py
bc50e523c208660fca1d44cd640de564e4046d434ba9c53d233f333fa8e35025  source/iterative/inversion.py
d76fcca13dec137e795d402a9b49ea7c01df8eacf5ed13136e0fe712e45bcd44  source/iterative/spectra.py
34e6ca4914207ab9da3596efd077fac34c788176263590e4159faed9d6c3cbd2  source/kinflation.py
693fee9368fb1a063dcf8f88e8be8b9167129914ea16cc666b522c6f42b050b8  source/latex_output.py
c6a3303bf62b5edc0980bed494235bb1594f7ebe5fd61e117efefd3db447bf67  source/mode_coefficients.py
05b1aa9da538813274b75a573d13b582e482c54ec7aba5f008dfc1a890677d36  source/orders.py
2c94e846155862778a2bf1bbed3a869a479ed0815a5a5b3e0fe4d610858b8326  source/presentation.py
00c9bd6accd6b0749f8f704a38af54ab30681e033965f7aefc377989ae4a4486  source/results.py
8fa2a012c5527f9431d3fe7e4d2a84ca30cb74237cdaee87cf4574b0efdc272a  source/series.py
66d11ab803fc7696f78f9a379121922402dfbd1ada748463d2da25e9703b4522  source/spectra.py
5f96456d959d2e746c0b6b7ecbf36c55397af36d9c752b5e8f3880208ad968be  source/terminal.py
58444ecaa045a66672d566b89ebf32cc436d83a106bbe54dbca227114654dd21  validation/canonical_SRI/LO/large_n_validation/run.py
56a97f0d043d2c0935a5aff9c769f85e3118cc04a26c5f1d1b9f9ff44df929e2  validation/canonical_SRI/LO/large_n_validation/archive/High_order_test.py
7e0e29ad3b6df130a206b181c5e02828bd86ac2bfedd5c03f5c617668612ccb0  validation/canonical_SRI/LO/large_n_validation/archive/Jet_test.py
c9b887209c062acaa7e013c7b411a88de99b852368096d1cc7853a5c88418343  validation/canonical_SRI/LO/large_n_validation/archive/Symbolic_test.py
05363dc7ad126ad2c0e7dc6ffe8ccf6325cf1a8ba62d817db2b1c1d566a1670e  validation/canonical_SRI/LO/large_n_validation/archive/projective_consistency_hierarchy 12.29.19.py
df63c0b17736cf6d3d963a065ae779465ccda38546a75304dc4b6b00535a62c8  validation/canonical_SRI/LO/large_n_validation/archive/projective_consistency_hierarchy_tests.py
```
