---
layout: post
title:  "Sequence alignment"
date:   2016-11-03 20:00
categories: jekyll
permalink: /2016/gotoh1
---

Sequence alignment arguably is one of the classic and fundamental problems in bioinformatics. It is not just that it has obvious applications (homology modelling, gene finding, evolutionary distances, ...) and [not so obvious applications](http://bioinformatics.oxfordjournals.org/content/26/18/i426.long), but the efficient algorithms for it (e.g. Needleman-Wunsch, Smith-Waterman, Gotoh, just to mention the most prominent ones) follow an algorithmic paradigm that can be applied to solve various optimization problems (to be precise, the ones following Richard Bellman's [Principle of Optimality](https://en.wikipedia.org/wiki/Bellman_equation#Bellman.27s_Principle_of_Optimality)). On top of that, one can learn a great deal about programming when implementing these algorithms. 

IMHO, every bioinformatics student (and probably everyone interested in doing computational work in biology) should implement these algorithms **correctly** and **efficiently**. To be fair, it is not exactly or not only my own opinion, but it is one of the main tasks in the practical course "Programmierpraktikum" that is mandatory for all bioinformatics students at LMU Munich. So, if you happen to encounter a bioinformatician that graduated in Munich, you can be quite sure that he or she implemented it. Of course, this also applies to myself. On top of that, since I started my PhD I have been responsible for the Sequence Alignment lecture that has been held for the students in the beginning of the practical course. So I guess this explains my enthusiasm for sequence alignment a bit...

## What is sequence alignment

There are many ways how to introduce sequence alignment, this is how I introduced it in my lectures (and I would also like to point you to the excellent script of Volker Heun's [Algorithmic bioinformatics I & II](https://www.bio.ifi.lmu.de/mitarbeiter/volker-heun/notes/index.html), were a different approach is used):

The most important biomolecules (DNA, RNA, proteins) have a primary structure that can be represented as a sequence of characters (in programming jargon also called a *string*). Every programmer knows how to answer the following question in his favourite programming language (I really, really hope that I am right about that!): *Are two strings equal*.

In the practical course, the students had to use Perl and Java (you might guess, all these Python guys were really, really angry about that...) so I also put the Perl and Java solutions onto the corresponding slide (using the eq operator in Perl and the equals method in Java) just to remind also the most slow-witted student not to use the == operator for that (yes I mean you, you laptop-open-but-not-making-notes-instead-checking-facebook-guy in the last row)... I am such a good teacher.

This of course is far away from what we actually want to ask in biology: *Are two strings similar*. Obviously, this is a much more difficult question, and does not belong to the standard library of any programming language (at least not that I know of). The difficulty starts with what is actually meant by *similarity* here. For human intuition, it is clear that WTHGQACVELSIW and WTHGRACVELSIW are quite similar (just an R for a Q in the second string). The strings WTHGQACVELSIW and WTHGACVELSIW are quite similar as well (now the Q is missing). But what about WTHGQACVELSIW and WTHAVSLW? What about the strings IIIIIIIIII and LLLLLLLLLL?

The most ambitious students even raised hands to answer the last question, of course correctly (naturally, as they are suitable for a university lecture in contrast to the laptop guy from the last row).

The lesson here is that *similarity* can mean different things, and it depends on what you are actually working on: An engineer at Google or Facebook trying to automatically find the advertisement you will fall for will use a different notion of similarity than a biologist. For a biologist, similarity is either related to the physico-chemical properties of the biomolecules, or their evolutionary history (naturally, these two notions are not independent of each other). 

Either way, the term *similarity* has to be defined formally. Usually, the notion of similarity comes down to a single number (a score), the higher this score for two sequences, the more similar they are.

## Alignments

For biological sequences, similarity is often defined via an alignment:

    WTHGQACVELSIW
    |||.     .|.|
    WTHA-----VSLW

Aligned characters (the ones with `|` or . in between) are called matches and mismatches, respectively (an evolutionary biologist may also speak of conserved positions and substitutions), the - are called gaps (or insertions/deletions).

Formally, an alignment for the sequences $$s_1...s_n$$ and $$t_1...t_m$$ is a sequence of integer pairs $$(x_i,y_i)$$ for $$i \in 1..k$$, such that 

1. $$x_i<x_{i+1}$$ and $$y_i<y_{i+1}$$
2. $$1 \leq x_i \leq n$$ and $$1 \leq y_i \leq m$$

The tuples denote all the aligned character positions. The formal alignment of the example above would therefore be

    (1,1),(2,2),(3,3),(4,4),(10,5),(11,6),(12,7),(13,8)

Intuitively, the `|` may not cross each other in the visual alignment, and a gap is never "aligned" with another gap.

## Similarity w.r.t. an alignment

The similarity score for a given alignment can now be broken down into a sum of similarity scores for each aligned position plus penalties (i.e. negative score components) for gaps. In the most basic case, matches are rewarded with some positive score (say, 3), mismatches are penalized with a negative score (say, -4). This is most often used when aligning DNA or RNA sequences. For proteins, a so-called [substitution matrix](https://en.wikipedia.org/wiki/Substitution_matrix) is used (derived from evolutionary considerations and data), where one can simply look up the similarity score for each pair of amino acids.

The most basic way to penalize gaps is to add a fixed gap penalty for each unaligned character (let's call that *fixed gap costs* for the moment, but we will return to that in the [next post](/2016/gotoh2)).

It is straight-forward to compute the score of a given alignment (and you definitely should implement that for testing your alignment algorithm implementation, see below). For the example above, it is just the sum

    11+5+8+0-2-2-2-2+1+4+2+11=31

when using the [BLOSUM62](https://de.wikipedia.org/wiki/BLOSUM) substitution matrix and gap costs of -2.

Now, let us look back what we actually wanted: We wanted to compute the similarity of two sequences. What we can do now is to compute the similarity of a specific alignment of two sequences. Does this help with our original problem? 

Well, of course: If we had a biological meaningful alignment of the two sequences in question, its similarity score would nicely reflect what we actually mean by similarity of the sequences. A meaningful alignment could for instance correspond to the structural superposition of two proteins or their evolutionary history.

But how do we get such a meaningful alignment? (And this is one of the central ideas behind sequence alignment, so I'll put it in bold face) **When you use a meaningful substitution matrix and gap costs, an alignment with maximal score among all possible alignments of two sequences is very likely also meaningful.**

So, (and I think this is not appreciated enough, so I stress this point here) the scoring system (substutution matrix and gap costs) is used for two distinct things in sequence alignment: to identify a meaningful alignment *and* to compute the alignment similarity. In fact, there may be applications, were different scoring systems should be used for these two things.

Now we finally can define our sequence alignment problem: Given two sequences, a substitution matrix and a gap penalty, compute the optimal alignment of the two sequences and its score.

## Naive algorithm

Now the supercomputing guy could shout out: I have a bazillion cores in my cluster, I just try every possible alignment of the two, compute their scores and keep the optimal one. Apparently, this guy is already bored by this blog article and wants to solve this problem by what he has already learned so far.

You might guess that this is not a practical approach. In fact, there is an astronomic number of 386 733 690 827 821 609 alignments for two sequences of length 20 and 30 (which would be ridiculously small proteins).

*Exercise: Give a formula to compute the number of alignments for given sequences of lengths m and n. Implement it and check the above mentioned number. What do I mean by the following statement: Here I used a slightly different definition of alignment.*

## Dynamic programming

As mentioned above, sequence alignment can be solved way more efficiently using dynamic programming (DP). The general principle of DP is as follows:

1. Break down your optimization problem into smaller subproblems.
2. Solve the smallest subproblems first and tabulate intermediate results.
3. Use tabulated results to solve greater subproblems until you end up with the solution of the whole problem.

Of course, this procedure is not applicable to each and every problem. It can only be used if one can prove that this indeed yields the optimal solution.

An important concept in DP is *backtracking*. Usually, intermediate results correspond to optimal scores for the subproblems. So the result of DP usually also is only the optimal score itself. However, by inspecting which intermediate results have been used to construct the final score, and inspecting which intermediate results have been used to construct these intermediate results and so on, one often can reconstruct the actual optimal solution (whatever that might be).

The basic idea of DP for sequence alignment is to solve sequence alignment for prefixes first (that would be our subproblems), tabulate the optimal alignment scores for prefixes and use them to quickly compute scores for longer prefixes, until the score for the whole sequences is computed. By backtracking the prefix alignments, the actual optimal alignment can be recovered.

For a sequence $$s_1...s_n$$, we denote its $$n+1$$ prefixes by $$s^0...s^n$$, where $$s^i=s_1...s_i$$.

Thus, our DP table for input sequences $$s_1...s_n$$ and $$t_1...t_m$$ will be a matrix with entries $$A_{i,j}$$ with $$0 \leq i \leq n$$ and $$0 \leq j \leq m$$ corresponding to the score of the optimal alignment of $$s^i$$ and $$t^j$$. Note that the matrix has size $$(n+1)$$ x $$(m+1)$$ corresponding to all combinations of the $$n+1$$ and $$m+1$$ prefixes, respectively.

Again (because this is so important): $$A_{i,j}$$ is the score of the optimal alignment of $$s^i$$ and $$t^j$$.

The basic cases are straight-forward: There is only one valid alignment of $$s^i$$ and $$t^0$$:

| $$s_1$$| $$s_2$$| $$s_3$$| ... | $$s_i$$ |
| - | - | - | ... | - |

and its score is $$A_{i,0} = g \cdot i$$ for a fixed gap cost g. Similarly $$A_{0,j} = g \cdot j$$. $A_{0,0}$ is not defined, but it is convenient to set it to zero (convenient means that this does not produce wrong results but prevents from introducing a special case for prefix alignments of length 1).

Now comes the most difficult part of this article: What about longer prefixes (and how do the tabulated values help here?):

From our alignment definition we know that the last column of the optimal alignment of $$s^i$$ and $$t^j$$ either contains (a) $$s_i$$ and $$t_j$$, or (b) $$s_i$$ and -, or (c) - and $$t_j$$. If (a) is the case:


| .............. | $$s_i$$ |
| .............. | $$t_j$$ |

Then, the alignment part (let's call it P) without the last column is an optimal alignment of $$s^{i-1}$$ and $$t^{j-1}$$.

**Prove**: First, we realize that P indeed is an alignent of $$s^{i-1}$$ and $$t^{j-1}$$. Assume that it is not optimal. Then there is an alignment Q of $$s^{i-1}$$ and $$t^{j-1}$$ with a better score than P. Replacing P by Q in the alignment above would result in a better score for $$s^{i}$$ and $$t^{j}$$, contradicting our assumption.

In other words: The dotted part of the alignment above also must be optimal (for shorter prefixes), because otherwise, a better alignment would increase the alignment score of the above alignment, which is impossible as it is optimal by assumption. Similar proves can be given for cases (b) and (c).

Thus, we can compute the scores for cases (a) to (c) and take the best one:

$$ A_{i,j} = \max\{ A_{i-1,j-1} + S(s_i,t_j), A_{i-1,j}+g, A_{i,j-1}+g\}$$

In tabular form, it becomes obvious how this algorithm proceedes:

| A | | $$t_1$$ | ... | $$t_{j-1}$$ | $$t_{j}$$ | ... |  $$t_m$$ |
| |
| $$s_1$$ |
| ... |
| $$s_{i-1}$$ ||||$$A_{j-1,j-1}$$|$$A_{i-1,j}$$
| $$s_i$$ |||| $$A_{i,j-1}$$ | $$A_{i,j}$$
| ... |
| $$s_n$$ |

An entry in the table can be filled when its top and left neighbors are already computed. The top-most and left-most row can be filled using the basic cases as mentioned above, so the algorithm has to start in the top left corner ($$A_{1,1}$$) and proceed, for instance row-by-row, until it ends up computing the bottom-right entry ($$A_{n,m}$$).

Obviously, the longest prefixes $$s^n$$ and $$t^m$$ are the whole strings $$s$$ and $$t$$, so $$A_{n,m}$$ is the score of the optimal alignment of $$s$$ and $$t$$. If you are additionally interested in the optimal alignment itself, it can be reconstructed using backtracking as mentioned above: Starting from the bottom-right entry (the one containing the optimal alignment score), proceed by always asking the question: How was this value computed and continue in the corresponding field in the table. The three possibilities arose due to the three cases (a) to (c), so it is straight-forward to reconstruct the actual optimal alignment.

This concludes the so-called [Needleman-Wunsch algorithm](https://www.ncbi.nlm.nih.gov/pubmed/5420325). Its runtime and space requirement obviously is in $$O(nm)$$.

## Variants

The type of alignment found by the Needleman-Wunsch algorithm is called a global alignment: Every character that is not aligned is penalized. For many applications, this is too restrictive:

1. For homology modelling, one might be interested whether there is a sequence in a database of domains with 3d structure that is similar to potentially only a part of a query sequence.
2. If you do not trust the domains in the database, it may even be wise to be even less restrictive and identify parts in database sequences that are similar to a part of the query sequence.
3. In many Next-Generation-Sequencing applications, the actual DNA/RNA fragment length is shorter than the read length (e.g. miRNA-seq, Ribosome profiling, PAR-CLIP). This makes it necessary to trim adapter sequences. Due to sequencing errors mismatches must be considered.

These three applications have in common that gaps at the beginning or end of one or both of the sequence should not be penalized: For 1. leading and trailing gaps in the database sequences should not be penalized, for 2. no leading or trainling gap should be penalized and for 3., a leading gap in the adapter sequence and a trailing gap in the read sequence should be free.

Can we reuse our good old Needleman-Wunsch algorithm for these special cases? Luckily, we can: For trailing gaps, just assume that the optimal alignment of two sequences $$s_1...s_n$$ and $$t_1...t_m$$ that does not penalize a trailing gap in $$s$$ is:

| .............. | $$s_n$$ | - | ... | - |
| .............. | $$t_j$$ | $$t_{j+1}$$ | ... | $$t_m$$ |

If the trailing gap is not penalized, the alignment score of $$s$$ and $$t$$ is equal to the alignment score of $$s$$ and $$t^j$$. Of course, this score is available in $$A_{j,m}$$. Thus, we only have to find the maximal score in the last row of $$A$$ to find the optimal alignment not penalizing a trailing gap in $$s$$. Similarly, we have to find the maximal score in the last column for optimal alignments not penalizing a trailing gap in $$t$$, or, the maximal value on both the last row and last column for the optimal alignment not penalizing a trailing gap in either $$s$$ or $$t$$.

But what if both should not be penalized? Such an optimal alignment looks like this:

| .............. | $$s_i$$ | - | ... | - |$$s_{i+1}$$ | ... | $$s_n$$ |
| .............. | $$t_j$$ | $$t_{j+1}$$ | ... | $$t_m$$ | - | ... | - |

The score of this optimal alignment (not penalizing the trailing gap in $$s$$ nor in $$t$$) is equal to the optimal global alignment score of $$s^i$$ and $$t^j$$. Thus, we only need to find the maximal value in the whole matrix $$A$$. In all cases, the optimal alignment itself can be reconstructed by backtracking starting from the corresponding entry in $$A$$.

But what about leading gaps that are supposed to be free? Let us assume that the optimal alignment not penalizing a leading gap in $$s$$ looks like this:

| - | ... | - | $$s_1$$ | .............. |
| $$t_1$$ | ... | $$t_{j-1}$$ | $$t_j$$ | .............. |

The alignment score in $$A$$ (whereever it might be due to special rules for trailing gaps) is too low by exactly $$\|(j-1) \cdot g\|$$. Thus, we have to redefine $$A$$: Previously, it contained the optimal alignment scores of prefixes, now it has to contain the alignment scores of prefixes *without penalizing gaps in $$s$$*. It is easy to see that this means to initialize $$A_{i,0}=0$$ for each $$i$$. For the analogue cases of not penalizing leading gaps in $$t$$ or not penalizing leading gaps in either $$s$$ or $$t$$, either the first row or both the first row and column must be initialized with 0.

In order not to penalize leading gaps in both sequences, $$A$$ must be redefined accordingly and computed in a slightly different manner: First, all basic cases must be initialized with 0. For all other entries, it is necessary to consider an additional case (d): The following alignment has score 0:

| - | ... | - | $$s_1$$ | ... | $$s_n$$ |
| $$t_1$$ | ... | $$t_m$$ | - | ...| - |

An alignment of this form also is a prefix of the full alignment that must not be penalized. Thus, the maximum must be computed as (the so-called *Waterman trick*)

$$ A_{i,j} = \max\{ A_{i-1,j-1} + S(s_i,t_j), A_{i-1,j}+g, A_{i,j-1}+g, 0\}$$


Importantly, these special cases for leading and trailing gaps can be combined, given full flexibility about which leading and trailing gaps are supposed to be penalized.

The special case where both leading and trailing gaps in both sequences are free is also called a local alignment and has first been brought forward by [Smith and Waterman](https://www.ncbi.nlm.nih.gov/pubmed/7265238). It runtime and space consumption still is $$O(nm)$$.	

## Overview

To compute sequence alignments for a given pair of sequences $$s_1...s_n$$ and $$t_1...t_m$$, a substitution matrix and gap costs, using dynamic programming the following steps must be done:

1. Create a matrix $$A$$ of size $$n+1$$ x $$m+1$$.
2. Initialize the first row and column
3. Use the recurrence formula above to systematically fill $$A$$
4. Report the similarity score
5. Use backtracking to reconstruct the alignment, if desired

The initialization and recurrence formula depends on the leading gap treatment, where to look for the similarity score depends on the trailing gap treatment.

