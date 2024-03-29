# Read in data functions

    Code
      read_contact_matrix(competition_path)
    Output
            B B.1 Ar Ar.1  M M.1  I I.1  E E.1  C C.1 Ai Ai.1  H H.1
      B     1   0  0    0  0   0  0   0  0   0  0   0 NA   NA NA  NA
      B.1   0   1  3    3  3   3 12  12  5   5  2   2 NA   NA NA  NA
      Ar   NA  NA  1    0  0   0  0   0  0   0  0   0  0    0  0   0
      Ar.1 NA  NA  0    1  2   2  3   3  4   4  2   2  3    3  1   1
      M    NA  NA NA   NA  2   0  1   0  0   0  0   1  0    0  0   0
      M.1  NA  NA NA   NA  0   2  3   4  2   2  3   4  3    3  3   3
      I    NA  NA NA   NA NA  NA 91   0  2  12  1   3  1    5  1   2
      I.1  NA  NA NA   NA NA  NA  0  91 10  24 33  37 21   27  8  11
      E    NA  NA NA   NA NA  NA NA  NA  3   0  1   0  0    1 NA  NA
      E.1  NA  NA NA   NA NA  NA NA  NA  0   3  4   5  7    8 NA  NA
      C    NA  NA NA   NA NA  NA NA  NA NA  NA 19   0  1    4  0   1
      C.1  NA  NA NA   NA NA  NA NA  NA NA  NA  0  19  3    8  5   6
      Ai   NA  NA NA   NA NA  NA NA  NA NA  NA NA  NA  4    0  0   4
      Ai.1 NA  NA NA   NA NA  NA NA  NA NA  NA NA  NA  0    4  1   5
      H    NA  NA NA   NA NA  NA NA  NA NA  NA NA  NA NA   NA  2   0
      H.1  NA  NA NA   NA NA  NA NA  NA NA  NA NA  NA NA   NA  0   2

---

    Code
      read_abundance(abundance_path, contact_matrix)
    Output
        species abundance
      1       B        27
      2      Ar        30
      3       M        23
      4       I       358
      5       E        62
      6       C       209
      7      Ai        95
      8       H        63

# Calculation interaction strength table

    Code
      interaction_strengths(contact_matrix, abundance, c(-0.1, -0.9, -0.2))
    Output
         Species_i Species_j Interactions Wins_i Wins_j Draws  F_ii        a_ii  F_ij
      1          B         B            1      0      0     1  -0.4 -0.01481481   0.0
      2          B        Ar            3      3      0     0   0.0  0.00000000  -0.3
      3          B         M            3      3      0     0   0.0  0.00000000  -0.3
      4          B         I           12     12      0     0   0.0  0.00000000  -1.2
      5          B         E            5      5      0     0   0.0  0.00000000  -0.5
      6          B         C            2      2      0     0   0.0  0.00000000  -0.2
      7         Ar        Ar            1      0      0     1  -0.4 -0.01333333   0.0
      8         Ar         M            2      2      0     0   0.0  0.00000000  -0.2
      9         Ar         I            3      3      0     0   0.0  0.00000000  -0.3
      10        Ar         E            4      4      0     0   0.0  0.00000000  -0.4
      11        Ar         C            2      2      0     0   0.0  0.00000000  -0.2
      12        Ar        Ai            3      3      0     0   0.0  0.00000000  -0.3
      13        Ar         H            1      1      0     0   0.0  0.00000000  -0.1
      14         M         M            2      0      0     2  -0.8 -0.03478261   0.0
      15         M         I            4      3      0     1   0.0  0.00000000  -0.5
      16         M         E            2      2      0     0   0.0  0.00000000  -0.2
      17         M         C            4      3      1     0   0.0  0.00000000  -1.2
      18         M        Ai            3      3      0     0   0.0  0.00000000  -0.3
      19         M         H            3      3      0     0   0.0  0.00000000  -0.3
      20         I         I           91      0      0    91 -36.4 -0.10167598   0.0
      21         I         E           24     10     12     2   0.0  0.00000000 -12.2
      22         I         C           37     33      3     1   0.0  0.00000000  -6.2
      23         I        Ai           27     21      5     1   0.0  0.00000000  -6.8
      24         I         H           11      8      2     1   0.0  0.00000000  -2.8
      25         E         E            3      0      0     3  -1.2 -0.01935484   0.0
      26         E         C            5      4      0     1   0.0  0.00000000  -0.6
      27         E        Ai            8      7      1     0   0.0  0.00000000  -1.6
      28         C         C           19      0      0    19  -7.6 -0.03636364   0.0
      29         C        Ai            8      3      4     1   0.0  0.00000000  -4.1
      30         C         H            6      5      1     0   0.0  0.00000000  -1.4
      31        Ai        Ai            4      0      0     4  -1.6 -0.01684211   0.0
      32        Ai         H            5      1      4     0   0.0  0.00000000  -3.7
      33         H         H            2      0      0     2  -0.8 -0.01269841   0.0
                  a_ij  F_ji        a_ji
      1   0.0000000000   0.0  0.00000000
      2  -0.0100000000  -2.7 -0.10000000
      3  -0.0130434783  -2.7 -0.10000000
      4  -0.0033519553 -10.8 -0.40000000
      5  -0.0080645161  -4.5 -0.16666667
      6  -0.0009569378  -1.8 -0.06666667
      7   0.0000000000   0.0  0.00000000
      8  -0.0086956522  -1.8 -0.06000000
      9  -0.0008379888  -2.7 -0.09000000
      10 -0.0064516129  -3.6 -0.12000000
      11 -0.0009569378  -1.8 -0.06000000
      12 -0.0031578947  -2.7 -0.09000000
      13 -0.0015873016  -0.9 -0.03000000
      14  0.0000000000   0.0  0.00000000
      15 -0.0013966480  -2.9 -0.12608696
      16 -0.0032258065  -1.8 -0.07826087
      17 -0.0057416268  -2.8 -0.12173913
      18 -0.0031578947  -2.7 -0.11739130
      19 -0.0047619048  -2.7 -0.11739130
      20  0.0000000000   0.0  0.00000000
      21 -0.1967741935 -10.6 -0.02960894
      22 -0.0296650718 -30.2 -0.08435754
      23 -0.0715789474 -19.6 -0.05474860
      24 -0.0444444444  -7.6 -0.02122905
      25  0.0000000000   0.0  0.00000000
      26 -0.0028708134  -3.8 -0.06129032
      27 -0.0168421053  -6.4 -0.10322581
      28  0.0000000000   0.0  0.00000000
      29 -0.0431578947  -3.3 -0.01578947
      30 -0.0222222222  -4.6 -0.02200957
      31  0.0000000000   0.0  0.00000000
      32 -0.0587301587  -1.3 -0.01368421
      33  0.0000000000   0.0  0.00000000

# Assemble a Jacobian matrix

    Code
      assemble_jacobian(interaction_table, abundance[, 1], "a_ij", "a_ji")
    Output
                  [,1]        [,2]         [,3]          [,4]         [,5]
      [1,] -0.01481481 -0.01000000 -0.013043478 -0.0033519553 -0.008064516
      [2,] -0.10000000 -0.01333333 -0.008695652 -0.0008379888 -0.006451613
      [3,] -0.10000000 -0.06000000 -0.034782609 -0.0013966480 -0.003225806
      [4,] -0.40000000 -0.09000000 -0.126086957 -0.1016759777 -0.196774194
      [5,] -0.16666667 -0.12000000 -0.078260870 -0.0296089385 -0.019354839
      [6,] -0.06666667 -0.06000000 -0.121739130 -0.0843575419 -0.061290323
      [7,]  0.00000000 -0.09000000 -0.117391304 -0.0547486034 -0.103225806
      [8,]  0.00000000 -0.03000000 -0.117391304 -0.0212290503  0.000000000
                    [,6]         [,7]         [,8]
      [1,] -0.0009569378  0.000000000  0.000000000
      [2,] -0.0009569378 -0.003157895 -0.001587302
      [3,] -0.0057416268 -0.003157895 -0.004761905
      [4,] -0.0296650718 -0.071578947 -0.044444444
      [5,] -0.0028708134 -0.016842105  0.000000000
      [6,] -0.0363636364 -0.043157895 -0.022222222
      [7,] -0.0157894737 -0.016842105 -0.058730159
      [8,] -0.0220095694 -0.013684211 -0.012698413

