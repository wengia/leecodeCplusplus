#ifndef Solution_h
#define Solution_h

#include<iostream>
#include <cmath>
#include <string>
#include <vector>
#include <stack>
#include <unordered_map>
#include <algorithm>
#include <queue>

using namespace std;

class Solution {
public:
	// Spiral Matrix II
    vector<vector<int> > generateMatrix(int n) {
        vector<vector<int>> res;
        vector<int> tmp;
        
        // Init the result;
        for(int i=0; i<n; i++) {
            tmp.clear();
            for(int j=0; j<n; j++)
                tmp.push_back(0);
            res.push_back(tmp);
        }
        
        // Start spiral matrix
		int imin = 0, imax = n-1; //row
        int jmin = 0, jmax = n-1; //col
		int count = 1;

		while(true) {
			for(int j=jmin; j<=jmax; j++)
				res[imin][j]=count++;
			if(++imin>imax) break;
			for(int i=imin; i<=imax; i++)
				res[i][jmax]=count++;
			if(jmin>--jmax) break;
			for(int j=jmax; j>=jmin; j--)
				res[imax][j]=count++;
			if(imin>--imax) break;
			for(int i=imax; i>=imin; i--)
				res[i][jmin]=count++;
			if(++jmin>jmax) break;
		}

        return res;
    }
};

#endif