# Errata

* **Page 227 (Eq. 10.96)** -- Missing minus sign in the first term on the right hand side (this error is also in the original paper Leondes et al., 1970). Thank you Filip Tronarp for spotting this. 

* **Page 265 (Alg. 12.10, Eq. 12.54)** -- The expression for P_0 does not hold in the general case (e.g., can be shown by the the Lyapunov equation FP + PF^T + LQL^T=0 not checking out). Thank you Vincent Adam for spotting this. After some further checking it has turned out that in fact P_0 is correct, but in general L and Q should be determined by L Q L^T = L_1 Q_1 L_1^T * P_2 + P_1 * L_2 Q_2 L_2^T, where * is the Kronecker product.

* **Page 260 (Eq. 12.34)** -- The power p should be D. Thank you Zheng Zhao for pointing this out.

* **Page 70 (Eq. 5.49)** -- The term m_u(v) should be m_u(t). Thank you Zheng Zhao for pointing this out.

* **Page 85 (Eq. 6.46)** -- Instead of \Sigma(\Delta), we should have \Sigma(\Delta t). Thank you Zheng Zhao for pointing this out.

* **Page 106 (Eq. 7.41)** -- Inconsistent notation, should read f^T(y, \tau) instead of f(y, \tau)^T (and the same for g). Thank you Zheng Zhao for pointing this out.

* **Page 227 (Eq. 10.98)** -- L(t) Q(t) L(t) P^-1(t) is missing a transpose on the second L(t). Thank you Zheng Zhao for pointing this out.
