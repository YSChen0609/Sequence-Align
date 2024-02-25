"""
Use Needleman-Wunsch algorithm and perform a pairwise sequence 
global alignment to find out their optimal alignment score, 
optimal alignment, and the corresponding sequence identity.

Note:
seq1 and seq2 have the tuple format: (seq_name, seq_string)
"""

import blosum
import pandas as pd
import numpy as np

class SeqAlign():

    def __init__(self, seq1, seq2, scoreFunc=blosum.BLOSUM(62)):
        self.seq1 = '*' + seq1[1]
        self.seq2 = '*' + seq2[1]
        self.scoreFunc = scoreFunc # default BLOSOM62
        
        # create a matrix for the Algo.
        self.m = len(self.seq1)
        self.n = len(self.seq2)
        self.M = [[None for _ in range(self.m)] for _ in range(self.n)]
        
        # create empty strings for optimal alignment
        self.align_seq1 = ''
        self.align_seq2 = ''
    
    def align(self):
        """
        Align query sequence to template sequence.
        Update the Matrix(self.M). 
        """
        
        # initialize the matrix
        self.M[0][0] = (self.scoreFunc['*']['*'], "END")
        
        ## deletions
        for i in range(1,self.m):
            self.M[i][0] = (self.M[i-1][0][0] + self.scoreFunc['*'][self.seq1[i]], "vertical")
            
        ## insertions
        for j in range(1,self.n):
            self.M[0][j] = (self.M[0][j-1][0] + self.scoreFunc['*'][self.seq2[j]], "horizontal")
        
        
        # Fill out the rest
        directions = ("diagonal", "vertical", "horizontal") # track back directions
        
        for i in range(1,self.m):
            for j in range(1,self.n):
                values = [
                    self.M[i-1][j-1][0] + self.scoreFunc[self.seq1[i]][self.seq2[j]], # substitution
                    self.M[i-1][j][0] + self.scoreFunc[self.seq1[i]]['*'],            # deletion
                    self.M[i][j-1][0] + self.scoreFunc['*'][self.seq2[j]]             # insertion
                ]
                
                # Find the best alignment score of the substrings and the direction
                max_val = max(values)
                direct  = directions[values.index(max_val)]
                
                self.M[i][j] = (max_val, direct)
        
    def getAlignment(self):
        """
        Update the self.align_seq and return the optimal alignment.
        """
        i = self.m - 1
        j = self.n - 1

        while self.M[i][j][1] != 'END': #(i>0 or j>0)
            if self.M[i][j][1] == 'diagonal':     # substitution
                self.align_seq1 = self.seq1[i] + self.align_seq1
                self.align_seq2 = self.seq2[j] + self.align_seq2
                i = max(i-1,0)
                j = max(j-1,0)
                
            elif self.M[i][j][1] == 'vertical':   # deletion
                self.align_seq1 = self.seq1[i] + self.align_seq1
                self.align_seq2 = '*' + self.align_seq2
                i = max(i-1,0)
                
            else:                                 # insertion
                self.align_seq1 = '*' + self.align_seq1
                self.align_seq2 = self.seq2[j] + self.align_seq2
                j = max(j-1,0)
        
        return (self.align_seq1, self.align_seq2)
    
    def getSeqIdentity(self):
        # Note that self.align() and self.getAlignment() have to be executed before calling this method.
        
        assert len(self.align_seq1) == len(self.align_seq2), "The aligned sequence should have the same length!"
        
        cnt = 0
        for (c1, c2) in zip(self.align_seq1, self.align_seq2):
            if c1 == c2 and c1 != '*':
                cnt += 1
        
        return (cnt/len(self.align_seq1))*100
        
        
    def printResult(self):
        # Note that self.align() and self.getAlignment() have to be executed before calling this method.
         
        print(f"The optimal alignment score: {self.M[-1][-1][0]}")
        print(f'The sequence identity is: {self.getSeqIdentity()}%\n')
        print(f"The optimal alignment:\n{self.align_seq1}\n\n{self.align_seq2}")
    
    def main(self):
        """
        Ger optimal alignment scores optimal alignment score,
        optimal alignment, and the corresponding sequence identity
        """
        self.align()
        self.getAlignment()
        self.printResult()