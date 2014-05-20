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
private:
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
		int m = word1.size(), n = word2.size();
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