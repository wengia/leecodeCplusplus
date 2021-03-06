#ifndef Solution_h
#define Solution_h

#include<iostream>
#include <cmath>
#include <string>
#include <vector>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <regex>
#include <assert.h>

using namespace std;

class Solution {
private:
	struct ListNode {
		int val;
		ListNode *next;
		ListNode(int x) : val(x), next(NULL) {}
	};

	struct TreeNode {
		int val;
		TreeNode *left;
		TreeNode *right;
		TreeNode(int x) : val(x), left(NULL), right(NULL) {}
	};
	TreeNode *root;

public:
	// For binary tree test!
	void createTree() {
		queue<TreeNode*> que;
		TreeNode *tmp;
		int x, ldata, rdata;
		int flag = -1;

		cout<<"\nCreate a new tree! Please imput the root ("<<flag<<" means null):";
		cin>>x;
		if (x==flag) {
			cout<<"Void tree\n";
			root = NULL;
			return;
		}
		root=new TreeNode(x);
		que.push(root);

		while (!que.empty()) {
			tmp=que.front();que.pop();
			cout<<"\nPlease imput the left and right child of node "<<tmp->val<<": ";
			cin>>ldata>>rdata;
			if(ldata!=flag) que.push(tmp->left=new TreeNode(ldata));
			if(rdata!=flag) que.push(tmp->right=new TreeNode(rdata));
		}
		cout<<"Successful creating a tree!\n";
	}

	// Balanced binary tree
	bool isBalanced() { return isBalanced(root); }

	bool isBalanced(TreeNode *root) {
		if(root==NULL) return true;
		
		int height = 0;
		return isBalanced(root, height);
    }

	bool isBalanced(TreeNode *root, int &height) {
		if(root==NULL) return true;

		int lh = 0, rh = 0;
		if(!isBalanced(root->left, lh)) return false;
		if(!isBalanced(root->right, rh)) return false;
		if(abs(lh-rh)>1) return false;

		height = 1 + (lh>rh? lh:rh);
		return true;
	}

	// Return the height of the tree
	int treeHeight() { return treeHeight(root); }

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

	//Add Binary
	string addBinary(string a, string b) {
        string res="";
		int count=0;

		for(int i=a.size()-1, j=b.size()-1; i>=0 || j>=0; i--, j--) {
			int add1 = (i<a.size())? (a[i]-'0'):0;
			int add2 = (j<b.size())? (b[j]-'0'):0;
			int sum = add1 + add2 + count;
			count = sum/2;

			if(sum==0 || sum==2) res.push_back('0');
			else res.push_back('1');
		}

		if(count==1) res.push_back('1');
		reverse(res.begin(), res.end());
		return res;
    }

	// Best Time to Buy and Sell Stock
	int maxProfit(vector<int> &prices) {
		int imin = 0, res = 0;

		for(int i=0; i<prices.size(); i++) {
			if(prices[i]<prices[imin])
				imin = i;
			res = max(res, prices[i]-prices[imin]);
		}

		return res;
    }

	// Best Time to Buy and Sell Stock II
	int maxProfitII(vector<int> &prices) {
		int imin = 0, tmp = 0, res = 0;

		for(int i=0; i<prices.size(); i++) {
			if(prices[i]<prices[imin])
				imin = i;
			if(tmp<prices[i]-prices[imin])
				tmp = prices[i]-prices[imin];
			else {
				res+=tmp;
				imin = i;
				tmp = 0;
			}
		}

		return res+tmp;
    }

	// Best Time to Buy and Sell Stock III
	int maxProfitIII(vector<int> &prices) {
		int n = prices.size(), imin = 0, imax = n-1, res = 0;
		if(n<2) return 0;

		int *l2r = new int[n];
		int *r2l = new int[n];
		
		l2r[imin] = 0;
		for(int i=1; i<n; i++) {
			if(prices[i]<prices[imin])
				imin = i;
			l2r[i] = max(l2r[i-1], prices[i]-prices[imin]);
		}

		r2l[imax] = 0;
		for(int i=n-2; i>=0; i--) {
			if(prices[i]>prices[imax])
				imax = i;
			r2l[i] = max(r2l[i+1], prices[imax]-prices[i]);
		}

		res = l2r[n-1];
		for(int i=0; i<n-1; i++)
			res = max(res, l2r[i]+r2l[i+1]);

		return res;
    }

	//Binary Tree Inorder Traversal
	vector<int> inorderTraversal() { return inorderTraversal(root); }
	vector<int> inorderTraversal(TreeNode *root) {
		struct Count{
			TreeNode *n;
			int count;
			Count(TreeNode *node, int c): n(node), count(c) {}
		};

        vector<int> res;
        stack<Count> que;
        
        if(root==NULL) return res;
        que.push(Count(root, 0));
        
		while(!que.empty()) {
			Count tmp = que.top(); que.pop();

			if(tmp.count==0) {
				tmp.count=1;
				que.push(tmp);
				if(tmp.n->left!=NULL) que.push(Count(tmp.n->left, 0));
			}
			else {
				res.push_back(tmp.n->val);
				if(tmp.n->right!=NULL) que.push(Count(tmp.n->right, 0));
			}
		}
        
        return res;
    }

	//Binary Tree Level Order Traversal 
	vector<vector<int> > levelOrder() { return levelOrder(root); }
	vector<vector<int> > levelOrder(TreeNode *root) {
		struct Node{
			TreeNode *tree;
			int level;
			Node(TreeNode *t, int l): tree(t), level(l) {}
		};

		vector<vector<int>> res;
		vector<int> floor;
		int level = 1;
		queue<Node> que;

		if(root==NULL) return res;
		que.push(Node(root, 1));

		while(!que.empty()) {
			Node tmp = que.front(); que.pop();

			if(tmp.level==level) floor.push_back(tmp.tree->val);
			else {
				level++;
				res.push_back(floor);
				floor.clear();
				floor.push_back(tmp.tree->val);
			}

			if(tmp.tree->left!=NULL) que.push(Node(tmp.tree->left, tmp.level+1));
			if(tmp.tree->right!=NULL) que.push(Node(tmp.tree->right, tmp.level+1));
		}
		res.push_back(floor);

		return res;
    }

	// Binary Tree Maximum Path Sum
	int maxPathSum() {
		int res = INT_MIN;
		maxPathSum(root, res);

		return res;
	}
	int maxPathSum(TreeNode *root, int &res) {
        if(root==NULL) return 0;

		int leftSum = maxPathSum(root->left, res);
		int rightSum = maxPathSum(root->right, res);
		int sum = max(root->val, root->val + max(leftSum, rightSum));
		res = max(res, sum);
		res = max(res, root->val + leftSum + rightSum);

		return sum;
    }

	//Binary Tree Zigzag Level Order Traversal
	vector<vector<int> > zigzagLevelOrder() {
        vector<vector<int>> res;
		zigzagLevelOrder(root, 0, res);

		for(int i=1; i<res.size(); i+=2)
			reverse(res[i].begin(), res[i].end());

		return res;
    }
	void zigzagLevelOrder(TreeNode *root, int level, vector<vector<int>> &res) {
        if(root==NULL) return;

		if(res.size()<=level)
			res.push_back(vector<int>());
		res[level].push_back(root->val);
		
		zigzagLevelOrder(root->left, level+1, res);
		zigzagLevelOrder(root->right, level+1, res);
    }

	// Candy
	int candy(vector<int> &ratings) {
		int *candy, res = 0;
		int n = ratings.size();

		candy = new int[n] ();
		for(int i=0; i<n; i++)
			candy[i] = 1;
		for(int i=1; i<n; i++)
			if(ratings[i]>ratings[i-1])
				candy[i] = candy[i-1]+1;
		for(int i=n-2; i>=0; i--)
			if(ratings[i]>ratings[i+1] && candy[i]<=candy[i+1])
				candy[i] = candy[i+1]+1;
		for(int i=0; i<n; i++)
			res+=candy[i];

		return res;
    }

	// Combinations
	vector<vector<int>> combine(int n, int k) {
        vector<vector<int>> res;
		combine(1, n, k, vector<int>(), res);
		return res;
    }

	void combine(int start, int n, int k, vector<int> tmp, vector<vector<int>> &res) {
		if(k==0) {
			res.push_back(tmp);
			return;
		}

		for(int i=start; i<=n; i++) {
			tmp.push_back(i);
			combine(i+1, n, k-1, tmp, res);
			tmp.pop_back();
		}
	}

	// Clone Graph
	struct UndirectedGraphNode {
		int label;
		vector<UndirectedGraphNode *> neighbors;
		UndirectedGraphNode(int x) : label(x) {};
	};
	UndirectedGraphNode *cloneGraph(UndirectedGraphNode *node) {
        unordered_map<UndirectedGraphNode *, UndirectedGraphNode *> mp;
		
		return cloneGraph(node, mp);
    }
    
    UndirectedGraphNode *cloneGraph(UndirectedGraphNode *node, unordered_map<UndirectedGraphNode *, UndirectedGraphNode *> &m ) {
        if(!node) return NULL;
		
		auto it = m.find(node);
		if(it!=m.end()) return it->second;

		UndirectedGraphNode *newNode = new UndirectedGraphNode(node->label);
		m[node] = newNode;
		for(int i=0; i<node->neighbors.size(); i++)
			newNode->neighbors.push_back(cloneGraph(node->neighbors[i], m ));
        
		return newNode;
    }

	// Construct Binary Tree from Inorder and Postorder Traversal
	void createTree(vector<int> inorder, vector<int> postorder) {
		root = buildTree(inorder, postorder);
	}

	TreeNode *buildTree(vector<int> &inorder, vector<int> &postorder) {
		int n = inorder.size();
		if(n==0) return NULL;
		TreeNode *tree = new TreeNode(postorder[n-1]);
		buildTree(inorder, postorder, 0, n-1, 0, n-1, tree);

		return tree;
    }

	void buildTree(const vector<int> &inorder, const vector<int> &postorder, 
		int in_start, int in_end, int post_start, int post_end, TreeNode *pre) {
		if(in_start>=in_end || post_start>=post_end) return;
		int in_pos;
		for(in_pos=in_start; in_pos<in_end; in_pos++)
			if(inorder[in_pos]==postorder[post_end]) break;
		int post_pos = post_start + in_pos - in_start - 1;

		// Add left child
		if(in_pos!=in_start) {
			pre->left = new TreeNode(postorder[post_pos]);
			buildTree(inorder, postorder, in_start, in_pos-1, post_start, post_pos, pre->left);
		}

		// Add right child
		if(in_pos!=in_end) {
			pre->right = new TreeNode(postorder[post_end-1]);
			buildTree(inorder, postorder, in_pos+1, in_end, post_pos+1, post_end-1, pre->right);
		}
	}

	// Construct Binary Tree from Inorder and Preorder Traversal
	void createTreePre(vector<int> preorder, vector<int> inorder) {
		root = buildTreePre(preorder, inorder);
	}
	TreeNode *buildTreePre(vector<int> &preorder, vector<int> &inorder) {
        int n = inorder.size();
		if(n==0) return NULL;
		TreeNode *tree = new TreeNode(preorder[0]);
		buildTreePre(preorder, inorder, 0, n-1, 0, n-1, tree);

		return tree;
    }

	void buildTreePre(vector<int> &preorder, vector<int> &inorder,
		int in_start, int in_end, int pre_start, int pre_end, TreeNode *pre) {
        if(in_start>=in_end || pre_start>=pre_end) return;

		int in_pos;
		for(in_pos=in_start; in_pos<in_end; in_pos++)
			if(inorder[in_pos]==preorder[pre_start]) break;
		int pre_pos = pre_end - (in_end - in_pos) + 1;

		// Add right child
		if(in_pos!=in_end) {
			pre->right = new TreeNode(preorder[pre_pos]);
			buildTreePre(preorder, inorder, in_pos+1, in_end, pre_pos, pre_end, pre->right);
		}

		// Add left child
		if(in_pos!=in_start) {
			pre->left = new TreeNode(preorder[pre_start+1]);
			buildTreePre(preorder, inorder, in_start, in_pos-1, pre_start+1, pre_pos-1, pre->left);
		}
    }

	//Convert Sorted Array to Binary Search Tree
	void createBST(vector<int> &num) {
		root = sortedArrayToBST(num);
	}
	TreeNode *sortedArrayToBST(vector<int> &num) {
        return sortedArrayToBST(num, 0, num.size()-1);
    }

	TreeNode *sortedArrayToBST(vector<int> &num, int first, int last) {
		if(first>last) return NULL;
		int mid = (first+last)/2;

		TreeNode *child = new TreeNode(num[mid]);
		child->left = sortedArrayToBST(num, first, mid-1);
		child->right = sortedArrayToBST(num, mid+1, last);

		return child;
	}

	// Copy List with Random Pointer
	struct RandomListNode {
		int label;
		RandomListNode *next, *random;
		RandomListNode(int x) : label(x), next(NULL), random(NULL) {}
	};
	RandomListNode *copyRandomList(RandomListNode *head) {
		if(!head) return NULL;
		unordered_map<RandomListNode*, RandomListNode*> m;
        
		RandomListNode dummy(0), *pre = &dummy, *tmp = head;
		while(tmp) {
			if(m.find(tmp)==m.end())
				m[tmp] = new RandomListNode(tmp->label);
			if(tmp->random && m.find(tmp->random)==m.end())
				m[tmp->random] = new RandomListNode(tmp->random->label);
			pre->next = m[tmp];
			pre = pre->next;
			if(tmp->random) pre->random = m[tmp->random];
		
			tmp = tmp->next;
		}

		return dummy.next;
    }

	// Distinct Subsequences
	int numDistinct_Answer(string S, string T) {
        int N = S.size(), M = T.size();
        int **dp;
		dp = new int*[M+1];
		for(int i=0; i<=M; i++){
			dp[i] = new int[N+1];
			dp[i][0] = 0;
		}
        dp[0][0] = 1;
        for (int j = 1; j <= N; ++j)
            dp[0][j] = 1;
        

        for (int i = 1; i <= M; ++i)
            for (int j = 1; j <= N; ++j)
                if (S[j-1] == T[i-1])
                    dp[i][j] = dp[i][j-1] + dp[i-1][j-1];
                else
                    dp[i][j] = dp[i][j-1];

        return dp[M][N];
    }

	int numDistinct(string S, string T) {
        int res = 0;
		numDistinct(S, T, 0, 0, res);
		return res;
    }

	void numDistinct(string S, string T, int start, int pos, int &res) {
		for(int i=start; i<S.size() && S.size()-i>=T.size()-pos; i++) {
			if(S[i]==T[pos]) {
				if(pos==T.size()-1) {
					res++;
					break;
				}
				numDistinct(S, T, i+1, pos+1, res);
			}

		}
	}

	//Edit Distance
	int minDistance(string word1, string word2) {
		const int m = word1.size(), n = word2.size();
		int **dp;

		dp = new int*[m+1];
		for(int i=0; i<=m; i++) {
			dp[i] = new int[n+1];
			dp[i][0] = i;
		}
		for(int j=0; j<=n; j++)
			dp[0][j]=j;

		for(int i=1; i<=m; i++)
			for(int j=1; j<=n; j++)
				if(word1[i-1]==word2[j-1])
					dp[i][j] = dp[i-1][j-1];
				else
					dp[i][j] = min(dp[i-1][j-1], min(dp[i-1][j], dp[i][j-1]))+1;

		return dp[m][n];
    }

	// Evaluate Reverse Polish Notation
	int evalRPN(vector<string> &tokens) {
        stack<int> str;
		int first, second;

		for(int i=0; i<tokens.size(); i++) {
			if(tokens[i].size()==1 && tokens[i]=="+") {
				second = str.top(); str.pop();
				first = str.top(); str.pop();
				str.push(first + second);
			}
			else if(tokens[i].size()==1 && tokens[i]=="-") {
				second = str.top(); str.pop();
				first = str.top(); str.pop();
				str.push(first - second);
			}
			else if(tokens[i].size()==1 && tokens[i]=="*") {
				second = str.top(); str.pop();
				first = str.top(); str.pop();
				str.push(first * second);
			}
			else if(tokens[i].size()==1 && tokens[i]=="/") {
				second = str.top(); str.pop();
				first = str.top(); str.pop();
				str.push(first / second);
			}
			else
				str.push(atoi(tokens[i].c_str()));
		}

		return str.top();
    }

	// Gas Station
	int canCompleteCircuit(vector<int> &gas, vector<int> &cost) {
		int total = gas.size();
        int res = 0, theMin = gas[0]-cost[0], sum = theMin;

		for(int i=1; i<total; i++) {
			sum += gas[i]-cost[i];
			if(sum<theMin) {
				theMin = sum;
				res = i;
			}
		}

		return sum>=0? (res+1)%total:-1;
    }

	// Grey Code
	vector<int> grayCode(int n) {
        vector<int> res(1<<n, 0);
		for(int i=0; i<1<<n; i++)
			res[i] = i ^ (i>>1);
		return res;
    }
	int grayToBinary(int gray) {
		for(int mask=gray>>1; mask>0; mask>>=1)
			gray ^= mask;
		return gray;
	}

	// Implement strStr()
	char *strStr(char *haystack, char *needle) {
        while(true) {
			char *tmp_h = haystack;
			char *tmp_n = needle;

			while(*tmp_n!='\0' && *tmp_h==*tmp_n) {
			    tmp_h++;
			    tmp_n++;
			}

			if(*tmp_n=='\0' ) return haystack;
			if(*tmp_h=='\0' ) return NULL;
			haystack++;
		}
    }

	// Insert Interval 
	struct Interval {
	   int start;
	   int end;
	   Interval() : start(0), end(0) {}
	   Interval(int s, int e) : start(s), end(e) {}
	};
	vector<Interval> insert(vector<Interval> &intervals, Interval newInterval) {
        bool insert = false;
		vector<Interval> res;

		for(auto it=intervals.begin(); it!=intervals.end(); it++) {
			if(insert || it->end<newInterval.start) {
				res.push_back(*it);
			} else if(newInterval.end<it->start) {
				res.push_back(newInterval);
				res.push_back(*it);
				insert = true;
			} else {
				newInterval.start = min(newInterval.start, it->start);
				newInterval.end = max(newInterval.end, it->end);
			}

		}
		if(!insert) res.push_back(newInterval);

		return res;
    }

	// Interleaving String
	bool isInterleave(string s1, string s2, string s3) {
		int size1=s1.size(), size2=s2.size();
		if(size1+size2!=s3.size()) return false;

		bool **dp;
		dp = new bool*[size1+1];
		dp[0] = new bool [size2 + 1];

		dp[0][0] = true;
		for (int i = 1; i <= size1; i++) {
			dp[i] = new bool[size2+1];
			dp[i][0] = dp[i - 1][0] && s1[i - 1] == s3[i - 1];
		}
		for (int j = 1; j <= size2; j++)
			dp[0][j] = dp[0][j-1] && s2[j-1] == s3[j-1];

		for (int i = 1; i <= size1; i++)
			for (int j = 1; j <= size2; j++)
				dp[i][j] = dp[i - 1][j] && s1[i - 1] == s3[i + j - 1] || 
							dp[i][j - 1] && s2[j - 1] == s3[i + j - 1];

		return dp[size1][size2];
    }

	// Largest Rectangle in Histogram
	int largestRectangleArea(const vector<int> &height) {
		//height.push_back(0);
		stack<int> st;
		int res = 0, i = 0, n = height.size();

		while (i < n) {
			if (st.empty() || height[i] >= height[st.top()])
				st.push(i++);
			else {
				int idx = st.top(); st.pop();
				int width = st.empty() ? i : i - st.top() - 1;
				res = max(res, width * height[idx]);
			}
		}

		return res;
	}

	// Length of Last Word
	int lengthOfLastWord(const char *s) {
		int count = 0, last = 0;

		while (*s != '\0') {
			if (*(s++) == ' ') {
				if (count != 0) last = count;
				count = 0;
			}
			else
				count++;
		}

		return count == 0 ? last : count;
	}

	// Longest Consecutive Sequence
	int longestConsecutive(vector<int> &num) {
		int res = 0, n = num.size();
		unordered_set<int> elements;

		for (int i = 0; i < n; i++)
			elements.insert(num[i]);
		for (int i = 0; i < n && !elements.empty(); i++) {
			int up = num[i], down = num[i];
			while (elements.find(up+1) != elements.end())
				elements.erase(++up);
			while (elements.find(down-1) != elements.end())
				elements.erase(--down);

			res = max(res, up-down+1);
		}

		return res;
	}

	// Longest Substring Without Repeating Characters
	int lengthOfLongestSubstring(string s) {
		int res = 0, size = s.size(), start = 0, end = 0;
		bool exist[256] = {false};

		while (end<size && start+res<size) {
			if (!exist[s[end]])
				exist[s[end++]] = true;
			else
				exist[s[start++]] = false;

			res = max(res, end - start);
		}

		return res;
	}

	// Maximal Rectangle, need to learn dp in this question
	int maximalRectangle(vector<vector<char> > &matrix) {
		if (matrix.size() == 0 || matrix[0].size() == 0) return 0;
		int m = matrix.size(), n = matrix[0].size();
		vector<int> height(n + 1, 0);
		int res = 0;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++)
				height[j] = (matrix[i][j] == '0') ? 0 : height[j] + 1;
			int tmp = largestRectangleArea(height);
			res = max(res, tmp);
		}
		
		return res;
	}

	// Maximum Subarray
	int maxSubArray(int A[], int n) {
		if (n == 0) return 0;
		int current, res;

		current = A[0];
		res = A[0];
		for (int i = 1; i < n; i++) {
			current = max(current + A[i], A[i]);
			res = max(res, current);
		}
		
		return res;
	}

	// Decode Ways
	int numDecodings(string s) {
		int n = s.size();
		if (n == 0 || s[0]=='0') return 0;
		int *dp = new int[n + 1];

		dp[0] = 1;
		dp[1] = 1;
		int ten, unit, isTen, isUnit;
		for (int i = 2; i <= n; i++) {
			unit = (s[i - 1] - '0');
			ten = (s[i - 2] - '0') * 10 + unit;
			isTen = (ten < 27 && ten > 9) ? 1 : 0;
			isUnit = (unit < 27 && unit > 0) ? 1 : 0;
			dp[i] = dp[i - 2] * isTen + dp[i - 1] * isUnit;
			if (dp[i] == 0) return 0;
		}

		return dp[n];
	}

	// Max Points on a Line
	struct Point {
		int x;
		int y;
		Point() : x(0), y(0) {}
		Point(int a, int b) : x(a), y(b) {}
	};

	void testMaxPoints() {
		vector<Point> points;
		points.push_back(Point(0, 0));
		points.push_back(Point(0, 0));
		points.push_back(Point(-1, -1));
		points.push_back(Point(2, 2));

		cout << maxPoints(points) << endl;
	}

	int maxPoints(vector<Point> &points) {
		unordered_map<double, int> line;
		int n = points.size();
		if (n == 0) return 0;

		int res = (n == 1) ? 1 : 0;
		
		for (int i = 0; i < n; i++) {
			line.clear();
			int repeat = 1;
			int local = 0;
			for (int j = i + 1; j < n; j++) {
				if (points[i].x == points[j].x && points[i].y == points[j].y) {
					repeat++;
					continue;
				}
				int divider = points[j].x - points[i].x;
				double slope = (divider == 0) ? numeric_limits<double>::infinity() : 1.0 * (points[j].y - points[i].y) / divider;
				line[slope]++;
				local = max(local, line[slope]);
			}
			res = max(res, local + repeat);
		}

		return res;
	}

	// Merge Intervals 
	static bool compare(Interval i1, Interval i2) { return i1.start < i2.start; }
	vector<Interval> merge(vector<Interval> &intervals) {
		vector<Interval> res;
		int size = intervals.size();
		if (size <= 1) return intervals;

		sort(intervals.begin(), intervals.end(), compare);
		int start = intervals[0].start, end = intervals[0].end;

		for (int i = 1; i < size; i++) {
			if (end < intervals[i].start) {
				res.push_back(Interval(start, end));
				start = intervals[i].start;
				end = intervals[i].end;
				continue;
			}
			if (end <= intervals[i].end)
				end = intervals[i].end;
		}
		res.push_back(Interval(start, end));

		return res;
	}

	// Merge Sorted Array
	void merge(int A[], int m, int B[], int n) {
		int i = m - 1;
		int j = n - 1;
		int k = m + n - 1;
		while (i >= 0 && j >= 0)
			if (A[i]>B[j])
				A[k--] = A[i--];
			else
				A[k--] = B[j--];
		while (j >= 0) A[k--] = B[j--];
	}

	// Merge Two Sorted Lists
	ListNode *mergeTwoLists(ListNode *l1, ListNode *l2) {
		ListNode dummy(0);
		ListNode *pre = &dummy;

		while (l1 && l2) {
			if (l1->val < l2->val) {
				pre->next = l1;
				l1 = l1->next;
			}
			else {
				pre->next = l2;
				l2 = l2->next;
			}
			pre = pre->next;
		}
		if (l1) pre->next = l1;
		if (l2) pre->next = l2;

		return dummy.next;
	}

	// Minimum Depth of Binary Tree
	int minDepth(TreeNode *root) {
		if (!root) return 0;
		if (!root->left && !root->right) return 1;

		int left = (!root->left)? INT_MAX : minDepth(root->left) + 1;
		int right = (!root->right) ? INT_MAX : minDepth(root->right) + 1;

		return min(left, right);
	}

	// Minimum Path Sum
	// Note: You can only move either down or right at any point in time.
	int minPathSum(vector<vector<int> > &grid) {
		if (grid.empty() || grid[0].empty()) return 0;
		int m = grid.size(), n = grid[0].size();
		int *dp = new int[n];

		dp[0] = grid[0][0];
		for (int i = 1; i < n; i++)
			dp[i] = dp[i - 1] + grid[0][i];

		for (int i = 1; i < m; i++) {
			dp[0] += grid[i][0];
			for (int j = 1; j < n; j++)
				dp[j] = min(dp[j - 1], dp[j]) + grid[i][j];
		}

		return dp[n - 1];
	}

	// Minimum Window Substring
	string minWindow(string S, string T) {
		int m = S.size(), n = T.size(), count = 0;
		int *need = new int[128] {0};
		int *find = new int[128] {0};
		int resstart = -1, resend = m;


		for (int i = 0; i < n; i++)
			need[T[i]] ++;

		for (int start = 0, end = 0; end < m; end++) {
			if (need[S[end]] == 0)
				continue;
			if (find[S[end]] < need[S[end]])
				count++;
			find[S[end]]++;
			if (count < n) continue;

			// find "start"
			for (; start < end; start++) {
				if (need[S[start]] == 0) continue;
				if (find[S[start]] == need[S[start]]) break;
				find[S[start]]--;
			}

			// get res
			if (end - start < resend - resstart) {
				resend = end;
				resstart = start;
			}
		}

		return (resstart == -1) ? "" : S.substr(resstart, resend - resstart + 1);
	}

	// N Queens
	vector<vector<string>> solveNQueens(int n) {
		vector<vector<string>> res;
		vector<pair<int, int>> board;
		vector<string> current;

		solveNQueens(n, board, current, res);

		return res;
	}

	void solveNQueens(const int n, vector<pair<int, int>> &board, vector<string> &current, vector<vector<string>> &res) {
		int row = current.size();
		if (row == n) {
			res.push_back(current);
			return;
		}

		current.push_back(string(n, '.'));
		for (int j = 0; j < n; j++) {
			if (!isSafe(board, row, j))
				continue;
			current[row][j] = 'Q';
			board.push_back(pair<int, int>(row, j));
			solveNQueens(n, board, current, res);

			current[row][j] = '.';
			board.pop_back();
		}
		current.pop_back();
	}

	bool isSafe(const vector<pair<int, int>> &board, int x, int y) {
		int b = board.size();

		for (int i = 0; i < b; i++) {
			if (x == board[i].first || y == board[i].second || abs(board[i].first - x) == abs(board[i].second - y))
				return false;
		}

		return true;
	}

	// Palindrome Partitioning
	vector<vector<string>> partition(string s) {
		vector<string> palindrome;
		vector<vector<string>> res;

		partition(s, 0, palindrome, res);

		return res;
	}

	void partition(const string s, int start, vector<string> &palindrome, vector<vector<string>> &res) {
		int end = s.size();
		if (start == end) {
			res.push_back(palindrome);
			return;
		}

		for (int length = 1; length + start <= end; length++) {
			int ins = 0;
			for (; ins < length / 2; ins++) {
				if (s[start + ins] != s[start + length - 1 - ins])
					break;
			}
			if (ins == length / 2) {
				palindrome.push_back(s.substr(start, length));
				partition(s, start + length, palindrome, res);
				palindrome.pop_back();
			}
		}
	}

	// Palindrome II
	int minCut(string s) {
		int n = s.size();
		if (n == 0 || n == 1) return 0;
		int *dp = new int[n];
		bool *isPalindrome = new bool[n];

		dp[0] = 0;
		for (int i = 1; i < n; i++) {
			dp[i] = dp[i - 1] + 1;
			isPalindrome[i] = true;
			for (int j = 0; j < i; j++) {
				isPalindrome[j] = (s[i] == s[j]) ? isPalindrome[j + 1] : false;

				if (isPalindrome[j])
					dp[i] = (j == 0) ? 0 : min(dp[i], dp[j - 1] + 1);
			}
		}

		return dp[n-1];
	}

	// Partition List
	ListNode *partition(ListNode *head, int x) {
		ListNode dummy(0), *pre = &dummy, *pos;
		dummy.next = head;

		// Find the place to insert
		while (pre->next != NULL && pre->next->val < x) pre = pre->next;

		// resort the list
		pos = pre;
		while (pos->next) {
			if (pos->next->val < x) {
				ListNode *tmp = pos->next;
				pos->next = pos->next->next;
				tmp->next = pre->next;
				pre->next = tmp;
				pre = pre->next;
			}
			else
				pos = pos->next;
		}

		return dummy.next;
	}

	// Pascal's Triangle
	vector<vector<int>> generate(int numRows) {
		vector<vector<int>> res;
		if (numRows == 0) return res;

		for (int i = 1; i <= numRows; i++) {
			vector<int> row(i, 1);
			for (int j = 1; j < i - 1; j++) {
				row[j] = res[i - 2][j - 1] + res[i - 2][j];
			}
			res.push_back(row);
		}

		return res;
	}

	// Pascal's Triangle II
	vector<int> getRow(int rowIndex) {
		vector<int> res(rowIndex + 1, 1);
		int before, now;

		for (int i = 2; i <= rowIndex; i++) {
			before = res[0];
			for (int j = 1; j < i / 2 + 1; j++) {
				now = res[j];
				res[j] = before + res[j];
				res[i - j] = res[j];
				before = now;
			}
		}

		return res;
	}

	// Permutation Sequence
	string getPermutation(int n, int k) {
		string res;
		int total = 1;
		for (int i = 1; i <= n; i++) {
			res.push_back(i + '0');
			total *= i;
		}

		//while (--k) nextPermutation(res);
		k--;
		while (n) {
			total /= n;
			int idx = k / total;
			k %= total;
			res.push_back(res[idx]);
			res.erase(idx, 1);
			n--;
		}

		return res;
	}

	void nextPermutation(string &num) { // Recursion Strategy
		int n = num.size(), pos;
		if (n <= 1) return;

		for (pos = n - 2; pos >= 0; pos--)
			if (num[pos] < num[pos + 1]) break;
			
		if (pos < 0) { // find max permutation
			reverse(num.begin(), num.end());
			return;
		}

		for (int i = n-1; i > pos; i--) {
			if (num[i] > num[pos]) {
				swap(num[i], num[pos]);
				sort(num.begin() + pos + 1, num.end());
				break;
			}
		}
	}

	// Permute
	vector<vector<int> > permute(vector<int> &num) {
		int n = num.size();
		vector<vector<int>> res;
		vector<int> current;
		vector<bool> avail(n, true);

		permute(num, avail, current, res);

		return res;
	}

	void permute(const vector<int> &num, vector<bool> &avail, vector<int> &cur, vector<vector<int>> &res) {
		int n = num.size();
		if (cur.size() == n) {
			res.push_back(cur);
			return;
		}

		for (int i = 0; i < n; i++) {
			if (!avail[i]) continue;
			cur.push_back(num[i]);
			avail[i] = false;
			permute(num, avail, cur, res);
			avail[i] = true;
			cur.pop_back();
		}
	}

	// Populating Next Right Pointers in Each Node
	struct TreeLinkNode {
		int val;
		TreeLinkNode *left, *right, *next;
		TreeLinkNode(int x) : val(x), left(NULL), right(NULL), next(NULL) {}
	};
	void connect(TreeLinkNode *root) {
		if (!root) return;
		queue<TreeLinkNode *> nodeQue;
		int level = 0, count = pow(2.0, level);

		nodeQue.push(root);
		while (!nodeQue.empty()) {
			TreeLinkNode *tmp = nodeQue.front(); nodeQue.pop();

			if (--count == 0) {
				tmp->next = NULL;
				level++;
				count = pow(2.0, level);
			}
			else {
				tmp->next = nodeQue.front();
			}

			if (tmp->left) nodeQue.push(tmp->left);
			if (tmp->right) nodeQue.push(tmp->right);
		}
	}

	// Populating Next Right Pointers in Each Node for any binary tree
	void connect_2(TreeLinkNode *root) {
		if (!root) return;
		TreeLinkNode *current = root, *last, *node;

		while (current) {
			last = NULL;
			node = current;
			while (node) {
				if (!last) current = node->left ? node->left : node->right;
				if (node->left || node->right) {
					if (last) last->next = node->left ? node->left : node->right;
					if (node->left && node->right) node->left->next = node->right;
					last = node->right ? node->right : node->left;
				}
				node = node->next;
			}
		}
	}

	// Recover Binary Search Tree
	void recoverTree(TreeNode *root) {
		if (!root || !root->left && !root->right) return;

		TreeNode *pre = NULL, *first = NULL, *second = NULL;
		recoverTree(pre, root, first, second);

		swap(first->val, second->val);
	}

	void recoverTree(TreeNode *&pre, TreeNode *cur, TreeNode *&first, TreeNode *&second) {
		if (!cur) return;

		if (cur->left) recoverTree(pre, cur->left, first, second);
		if (pre && pre->val > cur->val) {
			if (!first)
				first = pre;

			second = cur;
		}
		pre = cur;
		if (cur->right) recoverTree(pre, cur->right, first, second);
	}

	// Remove Duplicates from Sorted List 
	ListNode *deleteDuplicates(ListNode *head) {
		if (!head || !head->next) return head;
		ListNode *res = head, *pos = head->next;

		while (pos) {
			if (pos->val == head->val) {
				head->next = pos->next;
				delete pos;
				pos = head;
			}
			else
				head = head->next;
			pos = pos->next;
		}

		return res;
	}

	// Remove Duplicates from Sorted List II
	ListNode *deleteDuplicates_2(ListNode *head) {
		if (!head || !head->next) return head;
		bool dup;
		ListNode dummy(0), *pre = &dummy;
		dummy.next = head;

		while (head) {
			dup = false;
			while (head->next && head->val == head->next->val) {
				ListNode *tmp = head->next;
				head->next = head->next->next;
				delete tmp;
				dup = true;
			}
			if (dup) {
				pre->next = head->next;
				delete head;
				head = pre->next;
			}
			else {
				pre = pre->next;
				head = head->next;
			}
		}

		return dummy.next;
	}

	// Reorder List
	void reorderList(ListNode *head) {
		if (!head) return;
		ListNode *pos = head;
		stack<ListNode *> st;

		while (pos) {
			st.push(pos);
			pos = pos->next;
		}

		pos = head;
		while (pos != st.top() && pos->next != st.top()) {
			ListNode *rear = st.top(); st.pop();
			rear->next = pos->next;
			pos->next = rear;
			pos = pos->next->next;
		}
		st.top()->next = NULL;
	}

	// Restore IP Addresses
	vector<string> restoreIpAddresses(string s) {
		vector<string> res;
		if (s.size() < 4 || s.size() > 12) return res;
		string current;

		restoreIpAddresses(s, 0, 4, current, res);

		return res;
	}

	void restoreIpAddresses(const string s, int pos, int level, string &current, vector<string> &res) {
		int n = s.size();
		if (level == 0) {
			res.push_back(current);
			return;
		}

		int minIns = n - pos - (level - 1) * 3;
		int ins = minIns > 0 ? minIns : 1;
		for (;ins + pos <= n - level + 1 && ins <= 3; ins++) {
			string slot = s.substr(pos, ins);
			if (strcmp(slot.c_str(), "0")==0 || slot[0] != '0' && atoi(slot.c_str())<256) {
				if (level != 1) slot += '.';
				current.append(slot);
				restoreIpAddresses(s, pos + ins, level-1, current, res);
				current.erase(current.size() - slot.size(), slot.size());
			}
		}	
	}

	// Reverse Linked List II
	ListNode *reverseBetween(ListNode *head, int m, int n) {
		ListNode dummy(0), *pre = &dummy;
		dummy.next = head;
		
		n -= m;
		while (--m) pre = pre->next;
		head = pre->next;
		
		while (n--) {
			ListNode *move = head->next;
			head->next = move->next;
			move->next = pre->next;
			pre->next = move;
		}

		return dummy.next;
	}

	// Reverse Words in a String
	void reverseWords(string &s) {
		reverse(s.begin(), s.end());
		int start = s.size() - 1, end = s.size() - 1;

		while (start >= 0) {
			if (s[end] == ' ') {
				s.erase(end, 1);
				end--;
				start = end;
			}
			else if (s[start] == ' ') {
				reverse(s.begin() + start + 1, s.begin() + end + 1);
				end = --start;
			}
			else
				start--;
		}
		if (s[0] == ' ')
			s.erase(0, 1);
		else
			reverse(s.begin(), s.begin() + end + 1);
		
	}

	// Rotate List
	ListNode *rotateRight(ListNode *head, int k) {
		if (!head || k == 0) return head;
		int count = 1;
		ListNode *tail = head;

		while (tail->next) {
			tail = tail->next;
			count++;
		}
		k = k % count;
		if (k == 0) return head;

		ListNode *pre = head;
		count -= k;
		while (--count) pre = pre->next;

		tail->next = head;
		head = pre->next;
		pre->next = NULL;

		return head;
	}

	// Search a 2D Matrix
	bool searchMatrix(vector<vector<int> > &matrix, int target) {
		if (matrix.size() == 0 || matrix[0].size() == 0) return false;
		int m = matrix.size(), n = matrix[0].size(), i;

		for (i = 0; i < m; i++)
			if (target < matrix[i][0])
				break;
		if (i-- == 0) return false;
		
		for (int j = 0; j < n; j++)
			if (target == matrix[i][j])
				return true;

		return false;
	}

	// Set Matrix Zeroes
	void setZeroes(vector<vector<int> > &matrix) {
		if (matrix.size() == 0 || matrix[0].size() == 0) return;
		int m = matrix.size(), n = matrix[0].size();
		vector<int> row, col;

		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				if (matrix[i][j] == 0) {
					row.push_back(i);
					col.push_back(j);
				}

		for (int i = 0; i < row.size(); i++)
			for (int j = 0; j < n; j++)
				matrix[row[i]][j] = 0;

		for (int j = 0; j < col.size(); j++)
			for (int i = 0; i < m; i++)
				matrix[i][col[j]] = 0;
	}

	void setZeroes_2(vector<vector<int> > &matrix) {
		if (matrix.size() == 0 || matrix[0].size() == 0) return;
		int m = matrix.size(), n = matrix[0].size();
		int row = 0, col = 0;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (matrix[i][j] == 0) {
					row |= (1 << i);
					col |= (1 << j);
				}
			}
		}

		for (int i = 0; i < m; i++) {
			int tmp = row&(1 << i);
			for (int j = 0; tmp != 0 && j < n; j++)
				matrix[i][j] = 0;
		}
		for (int j = 0; j < n; j++) {
			int tmp = col&(1 << j);
			for (int i = 0; tmp != 0 && i < m; i++)
				matrix[i][j] = 0;
		}
	}

	// Search in Rotated Sorted Array II
	bool search(int A[], int n, int target) {
		int i = 0, j = n - 1;

		while (i < j && A[i] == A[j]) j--;
		while (i <= j) {
			int mid = (i + j) / 2;
			if (A[mid] == target)
				return true;
			if (A[i] <= A[mid]) {
				if (A[i] <= target && target < A[mid])
					j = mid;
				else
					i = mid + 1;
			}
			else {
				if (A[mid] < target && target <= A[j])
					i = mid;
				else
					j = mid - 1;
			}

		}

		return false;
	}

	// Simplify Path
	string simplifyPath(string path) {
		string res = "";
		path += "/";
		size_t pos, last = 0;

		pos = path.find_first_of("/");
		while (pos != string::npos) {
			string tmp = path.substr(last, pos - last);
			if (tmp == "..") {
				if (!res.empty())
					res.resize(res.find_last_of("/"));
			}
			else if (!tmp.empty() && tmp != ".") {
				res += ("/" + tmp);
			}
			last = pos + 1;
			pos = path.find_first_of("/", last);
		}

		return res.empty() ? "/" : res;
	}

	// Sort Colors
	void sortColors(int A[], int n) {
		int red = 0, blue = n-1, white = 0;
		
		while (white <= blue) {
			switch (A[white])
			{
			case 0:
				swap(A[white++], A[red++]);
				break;
			case 1:
				white++;
				break;
			case 2:
				swap(A[white], A[blue--]);
			}
		}
	}

	// Sqrt, Newton iteration
	int sqrt(int x) {
		if (x == 0) return x;
		double last = 0, current = 1;
		while (current != last) {
			last = current;
			current = (current + x / current) / 2;
		}

		return int(current);
	}

	// Subsets
	vector<vector<int> > subsets(vector<int> &S) {
		vector<vector<int>> res;
		vector<int> cur;

		sort(S.begin(), S.end());
		subsets_2(S, 0, cur, res);
		res.push_back(vector<int>());

		return res;
	}

	void subsets(const vector<int> &s, int start, vector<int> &cur, vector<vector<int>> &res) {
		for (int i = start; i < s.size(); i++) {
			cur.push_back(s[i]);
			res.push_back(cur);
			subsets(s, i + 1, cur, res);
			cur.pop_back();
		}
	}

	// Subsets II
	void subsets_2(const vector<int> &s, int start, vector<int> &cur, vector<vector<int>> &res) {
		for (int i = start; i < s.size(); i++) {
			if (i>start && s[i] == s[i - 1])
				continue;
			cur.push_back(s[i]);
			res.push_back(cur);
			subsets_2(s, i + 1, cur, res);
			cur.pop_back();
		}
	}

	// Sum Root to Leaf Numbers
	int sumNumbers(TreeNode *root) {
		int res = 0;
		sumNumbers(root, 0, res);
		return res;
	}

	void sumNumbers(TreeNode *root, int current, int &res) {
		if (!root) return;
		current = current * 10 + root->val;

		if (!root->left && !root->right) {
			res += current;
			return;
		}
		sumNumbers(root->left, current, res);
		sumNumbers(root->right, current, res);
	}

	// Symmetric Tree
	bool isSymmetric(TreeNode *root) {
		if (!root) return true;
		
		return isSymmetric(root->left, root->right);
	}

	bool isSymmetric(TreeNode *leftTree, TreeNode *rightTree) {
		if (!leftTree && !rightTree) return true;
		if (!leftTree || !rightTree) return false;

		return leftTree->val == rightTree->val &&
			isSymmetric(leftTree->left, rightTree->right) &&
			isSymmetric(leftTree->right, rightTree->left);
	}

	// Sudoku Solver
	void solveSudoku(vector<vector<char>> &board) {
		fillCells(board, 0, 0);
	}

	bool fillCells(vector<vector<char>> &board, int r, int c) {
		findNextBlank(board, r, c);
		if (r == 9 && c == 9) return true;
		for (int i = 1; i <= 9; i++) {
			if (!isAvail(board, r, c, '0' + i))
				continue;
			board[r][c] = '0' + i;
			if (fillCells(board, r, c + 1))
				return true;
			board[r][c] = '.';
		}

		return false;
	}

	void findNextBlank(const vector<vector<char>> &board, int &r, int &c) {
		int col;
		for (int i = r; i < 9; i++) {
			col = (i == r) ? c : 0;
			for (int j = col; j < 9; j++)
				if (board[i][j] == '.') {
				r = i;
				c = j;
				return;
				}
		}

		r = 9; c = 9;
	}

	bool isAvail(const vector<vector<char>> &board, const int r, const int c, const char num) {
		for (int i = 0; i < 9; i++)
			if (board[i][c] == num) return false;
		for (int j = 0; j < 9; j++)
			if (board[r][j] == num) return false;
		int cell_i = r / 3 * 3, cell_j = c / 3 * 3;
		for (int i = cell_i; i < cell_i + 3; i++)
			for (int j = cell_j; j < cell_j + 3; j++)
				if (board[i][j] == num) return false;
		return true;
	}

	// Text Justification
	vector<string> fullJustify(vector<string> &words, int L) {
		int pos = 0, letterCount = -1;
		vector<string> res;
		if ((words.empty() || words[0].empty()) && L==0) return res;

		for (int i = 0; i < words.size(); i++) {
			if (letterCount + 1 + words[i].size() > L) {
				// form a line
				res.push_back(formLine(words, L - letterCount + i - 1 - pos, pos, i - 1));

				// reset the variables
				letterCount = words[i].size();
				pos = i;
			}
			else {
				letterCount += (1 + words[i].size());
			}
		}
		res.push_back(words[pos]);
		int last = res.size() - 1;
		for (int i = pos + 1; i < words.size(); i++) {
			res[last] += (" " + words[i]);
		}
		res[last] += string(L - res[last].size(), ' ');

		return res;
	}

	string formLine(const vector<string> &words, int len, int start, int end) {
		string line;
		if (end == start) {
			line = words[start] + string(len, ' ');
			return line;
		}

		int odd = len % (end - start), ins = len / (end - start);
		line.append(words[start]);
		while (++start <= end) {
			line.append(odd-- > 0 ? 1 + ins : ins, ' ');
			line.append(words[start]);
		}

		return line;
	}

	// Triangle
	int minimumTotal(vector<vector<int> > &triangle) {
		int *dp, n = triangle.size();
		dp = new int[n];
		for (int row = 0; row < n; row++) {
			int before = dp[0];
			for (int j = 0; j <= row; j++) {
				if (j == 0 || j == row)
					dp[j] = before + triangle[row][j];
				else {
					int tmp = before;
					before = dp[j];
					dp[j] = min(dp[j], tmp) + triangle[row][j];
				}
			}
		}

		int res = dp[0];
		for (int i = 1; i < n; i++)
			res = min(res, dp[i]);
		return res;
	}

	// Unique Binary Search Trees
	int numTrees(int n) {
		if (n < 3) return n;
		int *dp = new int[n + 1];
		fill_n(dp, n + 1, 0);

		dp[0] = 1;
		dp[1] = 1;

		for (int i = 2; i <= n; i++) {
			for (int j = i - 1; j >= 0; j--)
				dp[i] += (dp[j] * dp[i - 1 - j]);
		}

		return dp[n];
	}

	// Unique Binary Search Trees II
	vector<TreeNode *> generateTrees(int n) {
		return generateTrees(1, n);
	}

	vector<TreeNode *> generateTrees(int value, int length) {
		vector<TreeNode *> children;
		if (length == 0) {
			children.push_back(NULL);
			return children;
		}
		for (int i = value; i < value + length; i++) {
			TreeNode *cur;
			vector<TreeNode *> left_children = generateTrees(value, i - value);
			vector<TreeNode *> right_children = generateTrees(i + 1, value + length - i - 1);
			
			for (int l = 0; l < left_children.size(); l++)
				for (int r = 0; r < right_children.size(); r++) {
					cur = new TreeNode(i);
					cur->left = copyTree(left_children[l]);
					cur->right = copyTree(right_children[r]);
					children.push_back(cur);
				}
		}

		return children;
	}

	TreeNode* copyTree(TreeNode *node) {
		if (!node) return NULL;
		TreeNode *c = new TreeNode(node->val);
		if (node->left) c->left = copyTree(node->left);
		if (node->right) c->right = copyTree(node->right);

		return c;
	}

	// Unique Path
	int uniquePaths(int m, int n) {
		if (m == 1 || n == 1) return 1;
		int path = m - 1 + n - 1;
		int choice = m > n ? (m - 1) : (n - 1);

		long long res = 1;
		for (int i = path; i > choice; i--)
			res *= i;
		for (int i = path - choice; i > 0; i--)
			res /= i;
		return res;
	}

	// Unique Path II
	int uniquePathsWithObstacles(vector<vector<int> > &obstacleGrid) {
		if (obstacleGrid.empty() || obstacleGrid[0].empty()) return 0;
		int m = obstacleGrid.size(), n = obstacleGrid[0].size();
		int **dp = new int*[m + 1];

		for (int i = 0; i <= m; i++) {
			dp[i] = new int[n + 1];
			dp[i][0] = 0;
		}
		for (int j = 1; j <= n; j++)
			dp[0][j] = 0;

		dp[0][1] = 1;
		for (int i = 1; i <= m; i++)
			for (int j = 1; j <= n; j++) {
				dp[i][j] = obstacleGrid[i - 1][j - 1] == 1 ? 0 : (dp[i - 1][j] + dp[i][j - 1]);
			}

		return dp[m][n];
	}

	// Valid Number
	bool isNumber(const char *s) {
		if (*s == '\0') return false;
		bool special = false;
		while (*s != '\0' && *s == ' ') s++;

		// start should not be 0 except for 0 and 0.xxxxxx
		if (*s == '\0' ||
			*s == 'e' || 
			*s == '.' && *(s + 1) == '\0' || 
			*s == '0' && *(s + 1) != '\0' && *(s + 1) != '.')
			return false;
		while (*s != '\0') {
			if (*s == '.' || *s == 'e') {
				if (special) return false;
				special = true;
			}
			else if (*s<'0' || *s>'9')
				break;
			s++;
		}
		while(*s!='\0') {
			if (*s != ' ') return false;
			s++;
		}

		return true;
	}

	// Valid Palindrome
	bool isPalindrome(string s) {
		if (s.empty()) return true;
		string s_cp;

		for (int i = 0; i < s.size(); i++)
			if (s[i] >= 'a' && s[i] <= 'z')
				s_cp.push_back(s[i]);
			else if (s[i] >= 'A' && s[i] <= 'Z')
				s_cp.push_back(s[i] - 'A' + 'a');

		int i = 0, j = s_cp.size() - 1;
		while (i < j) {
			if (s_cp[i++] != s_cp[j--]) return false;
		}

		return true;
	}

	// Validate Binary Search Tree
	bool isValidBST(TreeNode *root) {
		TreeNode *pre = NULL;
		
		return isValidBST(pre, root);
	}

	bool isValidBST(TreeNode *&pre, TreeNode *cur) {
		if (!cur) return true;
		if(!isValidBST(pre, cur->left)) return false;
		if (pre && pre->val >= cur->val)
			return false;
		pre = cur;
		return isValidBST(pre, cur->right);
	}

	// Word Break
	bool wordBreak(string s, unordered_set<string> &dict) {
		int n = s.size();
		bool *canBreak = new bool[n + 1];
		memset(canBreak, false, n + 1);

		canBreak[0] = true;
		for (int i = 1; i <= n; i++) {
			for (int j = i - 1; j >= 0; j--)
				if (canBreak[j] && dict.find(s.substr(j, i - j)) != dict.end()) {
					canBreak[i] = true;
					break;
				}
		}

		return canBreak[n];
	}

	// Word Break II
	vector<string> wordBreak_2(string s, unordered_set<string> &dict) {
		string current;
		vector<string> res;
		if (!wordBreak(s, dict)) return res;
		wordBreak_2(s, dict, 0, current, res);
		return res;
	}

	void wordBreak_2(const string &s, unordered_set<string> &dict, int start, string cur, vector<string> &res) {
		if (start == s.size()) {
			cur.pop_back();
			res.push_back(cur);
			return;
		}

		for (int i = start + 1; i <= s.size(); i++) {
			string tmp = s.substr(start, i - start);
			if (dict.find(tmp) != dict.end()) {
				wordBreak_2(s, dict, i, cur + tmp + ' ', res);
			}
		}
	}

	// Word Ladder
	int ladderLength(string start, string end, unordered_set<string> &dict) {
		queue<pair<int, string>> que;

		que.push(pair<int, string>(1, start));
		while (!que.empty()) {
			auto current = que.front(); que.pop();
			string word = current.second;
			for (size_t i = 0; i < word.size(); i++) {
				word = current.second;
				for (char j = 'a'; j <= 'z'; j++) {
					word[i] = j;
					if (word == end) return current.first + 1;
					if (dict.find(word) != dict.end()) {
						que.push(pair<int, string>(current.first + 1, word));
						dict.erase(word);
					}
				}
			}
		}

		return 0;
	}

	// Word Ladder II
	vector<vector<string>> findLadders(string start, string end, unordered_set<string> &dict) {
		vector<vector<string>> res;
		int length = INT_MAX;
		queue<pair<int, vector<string>>> que;

		que.push(pair<int, vector<string>>(1, vector<string>(1, start)));
		while (!que.empty()) {
			auto current = que.front(); que.pop();
			size_t size = current.second.size();
			string word = current.second[size-1];
			for (size_t i = 0; i < word.size(); i++) {
				word = current.second[size - 1];
				for (char j = 'a'; j <= 'z'; j++) {
					word[i] = j;
					vector<string> tmp(current.second);
					tmp.push_back(word);
					if (word == end) {
						if (current.first + 1>length) return res;
						length = current.first + 1;
						res.push_back(tmp);
					}
					if (dict.find(word) != dict.end()) {
						que.push(pair<int, vector<string>>(current.first + 1, tmp));
						dict.erase(word);
					}
				}
			}
		}

		return res;
	}

	// Word Search
	bool exist(vector<vector<char> > &board, string word) {
		if (board.empty() || board[0].empty() || word.empty()) return false;
		int m = board.size(), n = board[0].size();
		bool **avail = new bool*[m];

		for (int i = 0; i < m; i++) {
			avail[i] = new bool[n];
			for (int j = 0; j < n; j++)
				avail[i][j] = true;
		}

		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				if (board[i][j] == word[0] && isSubsequence(board, word, 0, i, j, avail))
					return true;

		return false;
	}

	bool isSubsequence(const vector<vector<char> > &board, const string &word, int word_pos, int r, int c, bool **avail) {
		if (!avail[r][c] || word[word_pos] != board[r][c]) return false;
		if (word_pos == word.size() - 1) return true;
		int m = board.size(), n = board[0].size();

		avail[r][c] = false;
		if (r>0 && isSubsequence(board, word, word_pos + 1, r - 1, c, avail)) return true;
		if (r<m - 1 && isSubsequence(board, word, word_pos + 1, r + 1, c, avail)) return true;
		if (c>0 && isSubsequence(board, word, word_pos + 1, r, c - 1, avail)) return true;
		if (c<n - 1 && isSubsequence(board, word, word_pos + 1, r, c + 1, avail)) return true;
		avail[r][c] = true;

		return false;
	}

	// Zigzag
	string convert(string s, int nRows) {
		if (nRows == 1) return s;
		string res;
		int ins = (nRows - 1) * 2, n = s.size(), i;

		for (int row = 0; row < nRows; row++) {
			i = 0;
			while (true) {
				if (row != 0 && row != nRows - 1) {
					int tmp = i * ins - row;
					if (tmp>0 && tmp < n) res.push_back(s[tmp]);
				}
				if (row + i*ins < n)
					res.push_back(s[row + (i++)*ins]);
				else
					break;
			}
		}

		return res;
	}

	// Wildcard Matching
	bool isMatch(const char *s, const char *p) {
		const char *sbackup = NULL, *pbackup;
		while (*s != '\0') {
			if (*s == *p || *p == '?') {
				s++;
				p++;
			} 
			else if (*p == '*') {
				while (*p == '*') p++;
				if (*p == '\0') return true;
				sbackup = s;
				pbackup = p;
			}
			else {
				if (sbackup == NULL) return false;
				s = ++sbackup;
				p = pbackup;
			}
		}
		while (*p == '*') p++;

		return *p == '\0';
	}

	// Integer to Roman
	string intToRoman(int num) {
		string res = "";
		
		if (num >= 1000) {
			int tmp = num / 1000;
			num -= tmp * 1000;
			while (tmp--) res.push_back('M');
		}
		if (num >= 100) res.append(intToOneLevelRoman(num, 100));
		if (num >= 10)  res.append(intToOneLevelRoman(num, 10));
		if (num >= 1)  res.append(intToOneLevelRoman(num, 1));

		return res;
	}

	string intToOneLevelRoman(int &num, int level) {
		string res = "", one, five, ten;
		if (level == 1) {
			one = "I"; five = "V"; ten = "X";
		}
		else if (level == 10) {
			one = "X"; five = "L"; ten = "C";
		}
		else if (level == 100) {
			one = "C"; five = "D"; ten = "M";
		}

		int tmp = num / level;
		num -= level * tmp;
		if (tmp == 9)
			res.append(one + ten);
		else if (tmp >= 5) {
			tmp -= 5;
			res.append(five);
			while (tmp--) res.append(one);
		}
		else if (tmp == 4)
			res.append(one + five);
		else
			while (tmp--) res.append(one);

		return res;
	}

	// Median of two sorted arrays
	double findMedianSortedArrays(int A[], int m, int B[], int n) {
		int total = m + n;
		if (total % 2 == 1)
			return findMedianSortedArrays(A, m, B, n, total / 2 + 1);
		else
			return (findMedianSortedArrays(A, m, B, n, total / 2) + findMedianSortedArrays(A, m, B, n, total / 2 + 1)) / 2;
	}

	double findMedianSortedArrays(int A[], int m, int B[], int n, int k) {
		if (m <= 0) return B[k - 1];
		if (n <= 0) return A[k - 1];
		if (k <= 1) return min(A[0], B[0]);

		if (m / 2 + n / 2 + 1 >= k) {
			if (A[m / 2] <= B[n / 2]) return findMedianSortedArrays(A, m, B, n / 2, k);
			else return findMedianSortedArrays(A, m / 2, B, n, k);
		}
		else {
			if (A[m / 2] <= B[n / 2]) return findMedianSortedArrays(A + m / 2 + 1, m - (m / 2 + 1), B, n, k - (m / 2 + 1));
			else return findMedianSortedArrays(A, m, B + n / 2 + 1, n - (n / 2 + 1), k - (n / 2 + 1));
		}
	}

private:
	// Return the height of one branch
	int treeHeight(TreeNode *tree) {
		if(tree==NULL) return 0;

		int l = treeHeight(tree->left);
		int r = treeHeight(tree->right);

		return 1 + ((l>r)? l:r);
	}
};

#endif