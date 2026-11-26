# Why the exact test works

Setting: $X_1, \dots, X_n$ are i.i.d. with a continuous distribution symmetric around $\theta$, and we test $H_0: \theta = 0$. Let $R_i$ be the rank of $|X_i|$ among $|X_1|, \dots, |X_n|$ and

$$T^+ = \sum_{i=1}^{n} R_i \, \mathbf{1}\{X_i > 0\}.$$

Continuity means there are no zeros and no ties with probability one.

## 1. Signs and magnitudes are independent

Under $H_0$, for any $t \ge 0$,

$$P(X > 0, |X| \le t) = P(0 < X \le t) = \tfrac12 P(-t \le X \le t) = P(X > 0)\,P(|X| \le t),$$

where the middle step uses symmetry. So $\mathrm{sign}(X_i)$ is independent of $|X_i|$, and by independence across observations the whole vector of signs is independent of the vector of magnitudes. Ranks and antiranks are functions of the magnitudes only, so they are independent of the signs too.

## 2. The null distribution is distribution free

Let $D_j$ be the antirank: the index of the observation with the $j$-th smallest magnitude. Define $W_j = \mathbf{1}\{X_{D_j} > 0\}$. Conditioning on $D = d$ and using step 1,

$$P(W = w \mid D = d) = P\big(\mathbf{1}\{X_{d_1} > 0\} = w_1, \dots, \mathbf{1}\{X_{d_n} > 0\} = w_n\big) = 2^{-n}.$$

The result does not depend on $d$, so $W_1, \dots, W_n$ are i.i.d. Bernoulli(1/2), and

$$T^+ = \sum_{j=1}^{n} j\, W_j.$$

The distribution of $T^+$ is therefore the same for every continuous symmetric $F$. Counting how many subsets of $\{1, \dots, n\}$ add up to $t$ gives $P(T^+ = t)$.

`signed_rank_pmf` computes it by convolution: if $p_{k}$ is the pmf using ranks $1..k$, then

$$p_k(t) = \tfrac12 p_{k-1}(t) + \tfrac12 p_{k-1}(t - k).$$

## 3. Symmetry and moments

The map $A \mapsto A^c$ is a bijection between subsets of $\{1, \dots, n\}$ with sum $t$ and subsets with sum $n(n+1)/2 - t$. Hence $T^+$ is symmetric around $n(n+1)/4$, which is also its mean. By independence,

$$\mathrm{Var}(T^+) = \sum_{j=1}^{n} \frac{j^2}{4} = \frac{n(n+1)(2n+1)}{24}.$$

## 4. Asymptotic normality

Write $T^+ - E[T^+] = \sum_j j\,(W_j - \tfrac12)$, a weighted sum of i.i.d. centred bounded variables with weights $a_j = j$. The Lindeberg condition holds because no single weight dominates:

$$\frac{\max_j |a_j|}{\sqrt{\sum_j a_j^2}} = \frac{n}{\sqrt{n(n+1)(2n+1)/6}} \longrightarrow 0.$$

Therefore

$$\frac{T^+ - n(n+1)/4}{\sqrt{n(n+1)(2n+1)/24}} \xrightarrow{d} N(0, 1).$$

## 5. Comparing discrete tests fairly

$T^+$ and the sign statistic are discrete, so a non-randomised test usually cannot reach a size of exactly $\alpha$. With $n = 12$ and $\alpha = 0.05$, the sign test can only achieve 0.019 and the signed-rank test 0.046. Comparing their raw power against a t test of size 0.05 penalises them for that gap. The power study uses the randomised version: reject if $T > k$, and with probability $\gamma$ if $T = k$, where

$$P_{H_0}(T > k) + \gamma\, P_{H_0}(T = k) = \alpha.$$

In the simulations we average the rejection probability itself, which has lower variance than drawing the extra coin flip.

## References

- Wilcoxon, F. (1945). Individual comparisons by ranking methods. *Biometrics Bulletin*, 1(6), 80-83.
- Hettmansperger, T. P. and McKean, J. W. (2011). *Robust Nonparametric Statistical Methods*, 2nd ed. CRC Press.
- Lehmann, E. L. (1999). *Elements of Large-Sample Theory*. Springer.
