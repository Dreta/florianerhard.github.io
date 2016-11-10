---
layout: post
title:  "The algorithm of Gotoh"
date:   2016-11-04 20:00
categories: jekyll
permalink: /2016/gotoh2
---

In the [last article](/2016/gotoh1), I introduced the sequence alignment problem. What I did not tell you is that it is practically useless for any real-world application. Oops, I am sorry about that!

But don't be upset, reading it was not a complete waste of your precious time: The algorithm of Gotoh I am going to present here is much easier to comprehend when you have a good knowledge about Needleman-Wunsch, and what you can do with Gotoh indeed is highly relevant in many applications.

## Why are *fixed gap costs* way too simple?

Remember our [example alignment](/2016/gotoh1#alignments)?

    WTHGQACVELSIW
    |||.     .|.|
    WTHA-----VSLW

Using the BLOSUM50 matrix and gap costs of $$g=-2$$, this alignment has a score of 40.9. Unfortunately, there is another alignment with a score of 45.9:

    WTHGQACVELSIW
    |||  | |  |.|
    WTH--A-V--SLW

Naturally, this alignment has a better score, as it only has a single mismatch as compared to three. In fact, it is the optimal alignment for these two sequences and this scoring system (BLOSUM50 and $$g=-2$$).

But remember what I [told you](/2016/gotoh1#similarity-wrt-an-alignment) about the importance of a meaningful alignment? 
Substitutions between glycine and alanine, leucine and valine or leucine and isoleucine occur frequently due to their biochemical similarities. Apart from that, there is only a single insertion or deletion (of five amino acids) in the first alignment. In the second, there are three insertions or deletions (of two, one and two amino acids, respectively). Judge for yourself, which alignment is more meaningful? Of course, noone can really know about the true evolutionary history of these two sequences (I am sure about that, because it was me who invented these two as an example for this blog), but in general, fewer gaps, even if they are long, should be preferred over many short gaps.

Unfortunately, using *fixed gap costs* does not allow us to model the length of gaps: Each gap by itself costs a fixed penalty, independent of other adjacent gaps. Thus, we have to adapt our scoring system accordingly: Instead of fixed gap costs, we now introduce a *gap cost function* $$g(n)$$, which maps the natural number $$n=1,2,...$$ to a gap penalty. 

Of note, the aforementioned fixed gap costs correspond to a linear gap cost function without intercept term:

$$ g(n) = -2 \cdot n$$

But let us now try a different class of functions for $$g$$, for instance the affine function $$g(n) = -10 - 2 \cdot n$$ and see what happens: The total gap cost for the first alignment is $$g(5) = -10 -2 \cdot 5 = -20$$. For the second alignment $$g(2)+g(1)+g(2)=-10-2 \cdot 2 -10-2 \cdot 1 -10-2 \cdot 2 = -40$$. Thus, the total score of the first alignment becomes 30.9, the total score of the second 15.9. In fact, the first alignment is optimal, and this is exactly what we wanted.

Thus, using an affine function $$g(n) = o + e \cdot n$$ can help to make the optimal alignments more meaningful. Of note, the parameters $$o$$ and $$e$$ are often called *gap open* and *gap extension* costs for obvious reasons. Be careful, when defined in a straight-forward way as the parameters of an affine function, for a gap of length 1, gap open **and** gap extension must be paid. This is a bit counterintuitive so many people use a different definition (such that a gap of length 1 only costs one gap open, of length 2 a gap open and a gap extension and so on). Of course, both definitions are equivalent, but when specifying parameters, always make sure that you know which definition has been used.

## Dynamic programming for sequence alignment with arbitrary gap cost functions

Can we still use our good old Needleman-Wunsch implementation to optimize alignments with any gap cost function? Unfortunatly, the answer is no, as the following example shows:

Remember our three cases for the dynamic programming algorithm with fixed gap costs (a.k.a. linear gap cost function)? Case (b) looked like that:

| .............. | $$s_i$$ |
| .............. | - |

With linear gap costs, if the optimal alignment of prefixes $$s^i$$ and $$t^j$$ looks like that, then the alignment part without the last column (the dotted part) also is an optimal alignment (for $$s^{i-1}$$ and $$t^j$). For other gap cost functions, this does not hold anymore:

Consider the optimal alignment of `AILW` and `AL` with $$g(n)=-10+2n$$:

    AILW
    AL--

As the similarity score of I and L is not much worse than the score for L and L, the affine gap cost function prefers one gap of length two over two gaps of length one. This alignment corresponds to case (b), however, the alignment

    AIL
    AL-

obviously is not optimal. Oh no! Do we have to start over and come up with something completely new? Luckily that's not the case: What leads to problems here is a long trailing gap in case (b) or (c). Therefore, we have to distinguish again three cases to derive how to compute $$A_{i,j}$$: (a) the optimal alignment of $$s^i$$ and $$t^j$$ ends with a match/mismatch, (b) the optimal alignment ends with a gap of length k in $$t$$, or (c) the optimal alignment ends with a gap of length k in $$t$$. Note that this is very similar to the derivation for linear gap cost functions, we only consider the trailing gap length and introduce an open variable k.

Again, in all three cases, there is a prefix alignment that is optimal. For (a), it is the alignment without the last column, which can be proven as for linear gap costs. For (b) (and, equivalently, for (c)), the optimal alignment has the form:

| .............. | x | $$s_{i-k+1}$$ | ... | $$s_i$$ |
| .............. | $$t_j$$ | - | ... | - |

x is either a gap or $$s_{i-k}$$. The alignment left of (and including) x and $$t_j$$ (let's call it P) is optimal under all alignments ending with $$t_j$$ (i.e. not having a gap in the last column in $$t$$).

**Prove**: Assume that it is not optimal. Then there is an alignment Q of $$s^{i-k}$$ and $$t^{j}$$ with a better score than P but still ending with $$t_j$$. Replacing P by Q in the alignment above would result in a better score for $$s^{i}$$ and $$t^{j}$$, contradicting our assumption.

Note that considering $$t_j$$ here as the alignment column before the trailing gaps in $$t$$ is crucial, as this ensures that there is no trailing gap in the shorter alignment (also see the counter-example above). By the character $$t_j$$, we ensure that the scores of the prefix alignment (either P or Q) and the trailing gap can be summed for a total alignment score.

Unfortunately, this makes things a bit more complicated, as we do not know k. We only know that, if (b) is the case, one of the scores in $$A_{0,j}$$ through $$A_{i-1,j}$$ plus the corresponding costs for the trailing gap leads to an optimal alignment for $$s^i$$ and $$t^j$$. Thus, to compute $$A_{i,j}$$ we now have to consider much more entries than before:

$$
\begin{align*} 
	A_{i,j}  & = \max\{ A_{i-1,j-1} + S(s_i,t_j), D_{i,j}, I_{i,j}\}\\
	I_{i,j}  & = \max_{1 \leq k \leq i}\{ A_{i-k,j}+g(k)\}\\
	D_{i,j}  & = \max_{1 \leq k \leq j}\{ A_{i,j-k}+g(k)\}
\end{align*} 
$$

For convenience, we factorized the maximum into the three cases (a), (b) and (c). In tabular form:

| A | | $$t_1$$ | ... | $$t_{j-1}$$ | $$t_{j}$$ | ... |  $$t_m$$ |
| | | | | | $$A_{0,j}$$
| $$s_1$$ |||||$$A_{1,j}$$
| ... |||||...
| $$s_{i-1}$$ ||||$$A_{j-1,j-1}$$|$$A_{i-1,j}$$
| $$s_i$$ |$$A_{i,0}$$|$$A_{i,1}$$|...| $$A_{i,j-1}$$ | $$A_{i,j}$$
| ... |
| $$s_n$$ |

I.e. we have to maximize over the whole column above the entry (and add, to each entry, the respective gap costs) and the whole row left if the entry. This increases the runtime to $$O(nm\max(n,m))$$ and is known as the algorithm of Waterman, Smith and Beyer.



## The algorithm of Gotoh for affine gap cost functions

Remember that in contrast to the linear gap cost function $$g(n) = g \cdot n$$, more general gap cost functions can produce more meaningful optimal alignments. However, these better alignments come with the cost of an often prohibitive cubic runtime.

However, can we do better when we restrict ourselves to the affine gap cost function $$g(n) = o + e \cdot n$$? Let's see:

We start by inserting that into the the definition of $I$ and factorize:

$$
\begin{align*} 
I_{i,j}  & = \max_{1 \leq k \leq i}\{ A_{i-k,j}+g(k)\}\\
& = \max_{1 \leq k \leq i}\{ A_{i-k,j}+o + e \cdot k \} \\
& = \max\{ A_{i-1,j}+o + e, \max_{2 \leq k \leq i}\{ A_{i-k,j}+o + e \cdot k \}\}\\
\end{align*} 
$$

Adjusting the index and factoring out an $$e$$ from the maximum, this becomes

$$
\begin{align*} 
&  = \max\{ A_{i-1,j}+o + e, \max_{1 \leq k \leq i-1}\{ A_{i-1-k,j}+o + e \cdot k \} +e \}\\
 & = \max\{ A_{i-1,j}+o + e, I_{i-1,j} +e \}
\end{align*} 
$$

What does that mean? If, in addition to $$A$$ we also tabulate $$D$$ and $$I$$, we can compute each $$I$$ and $$D$$ value in constant time! So, our overall runtime again is $$O(nm)$$, we only have to take care of three matrices instead of one during dynamic programming. I called them $$D$$ and $$I$$, as their scores correspond to the optimal alignments that end with a gap in $$s$$ (deletion) or in $$t$$ (insertion), respectively.

This is called the algorithm of [Gotoh](http://www.sciencedirect.com/science/article/pii/0022283682903989).

## Gotoh in a nutshell

First, all three matrices are initialized:

$$
\begin{align*} 
A_{0,k} &  = g(k)\\
A_{k,0} &  = g(k)\\
I_{0,k} &  = -\mbox{Inf}\\
D_{k,0} &  = -\mbox{Inf}
\end{align*} 
$$

For convenience, the undefined values $$I_{0,k}$$ and $$D_{k,0}$$ are set to infinity (They are undefined, as there is no alignment with an empty string $$s$$ that ends with gap in $$t$$): To compute $I_{1,1}$, for instance, the formula is $$I_{1,1} = \max\{ A_{0,1}+o + e, I_{0,1} +e \}$$. As $$I_{0,1} = -\mbox{Inf}$$, this undefined value will not be used and no special case must be introduced.

Second, fill the three matrices systematically with the recurrences:

$$
\begin{align*} 
	I_{i,j}  & = \max\{ A_{i-1,j}+o + e, I_{i-1,j} +e \}\\
	D_{i,j}  & = \max\{ A_{i,j-1}+o + e, D_{i,j-1} +e \}\\
	A_{i,j}  & = \max\{ A_{i-1,j-1} + S(s_i,t_j), D_{i,j}, I_{i,j}\}
\end{align*} 
$$

This can be done again e.g. row-by-row, always computing $$I$$, then $$D$$ and then $$A$$ and then proceeding to the next entry. The optimal alignment score can again be found in $$A_{n,m}$$.

Third, the optimal alignment itself can be computed using backtracking. Importantly, in a canonical implementation, not only the two coordinates must be tracked, but also in which matrix the current backtracking step is taking place. This always was the greatest cause for programming errors in our practical course, so we came up with a more clever backtracking scheme (see the [next](/2016/gotoh3) post for that; o yes, this is just a perfect cliffhanger...)

One important thing to consider is the following: What to do, if the maximum is not unique, so there are two (or more) possible ways for backtracking. The answer is clear: It does not matter, both alignments have the same optimal score!


## Additional remarks

Of course, all variants for special treatment of leading and trailing gaps (including local alignments) can also be computed with the algorithm of Gotoh (you just need to use the Waterman trick in $$A$$).

Of note, for concave gap cost functions, there is an [algorithm](http://link.springer.com/article/10.1007/BF02459948) with runtime $$O(nm \log(n+m))$$. 
