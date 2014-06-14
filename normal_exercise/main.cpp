// for exercise

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <stack>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include "solution.h"

using namespace std;

struct ListNode {
	int val;
    ListNode *next;
    ListNode(int x) : val(x), next(NULL) {}
};


void output( ListNode *head ) {
	ListNode *tmp=head;
	
	while( tmp!=NULL ) { cout<<tmp->val<<' '; tmp=tmp->next; }
}


vector<int> twoSum(vector<int> &numbers, int target) {
    int i, j;
    bool flag = false;
    vector<int> result;
        
    for ( i=0; i<numbers.size(); i++ ) {
        for ( j=i+1; j<numbers.size(); j++ ) {
            if ( numbers[i]+numbers[j]==target ) {
                result.push_back(i);
                result.push_back(j);
                flag = true;
                break;
            }
        }
            
        if (flag) break;
    }
        
    return result;
}

int longestNonRepeatString(string s){
	string tmp_1, tmp_2;
	int idx_1, idx_2, max=0, count, tmpmax;

	for (int i=0; i<s.length(); i++){
		tmpmax = 0;

		for (int j=0; j<i; j++){
			count = 0;
			if ( s[i]==s[j] ) {
				idx_1=i+1;
				idx_2=j+1;

				while ( idx_1<i && idx_2<s.length() && s[idx_1++]==s[idx_2++] );
				count = (j>idx_2-j)? j:(idx_2-j);
			}
			else
				tmpmax++;

			tmpmax = (tmpmax>count)? tmpmax:count;
		}

		max = ( max>tmpmax )? max:tmpmax;
	}

	return max;
}

string longestPalindrome(string s) {
    int res = 1, ins, pos=0;
        
    if ( s.length()<3 ) return s;
        
    for (int i=1; i<s.length(); i++ ) {
        ins = 1;

		//odd
        while ( i-ins>=0 && i+ins<s.length() && s[i-ins]==s[i+ins] ) ins++;
        if ( (ins-1)*2+1>res ) {
            res = (ins-1)*2+1;
            pos = i-ins+1;
        }

		//even
		ins = 1;
		while ( i-ins>=0 && i+ins-1<s.length() && s[i-ins]==s[i+ins-1] ) ins++;
        if ( (ins-1)*2>res ) {
            res = (ins-1)*2;
            pos = i-ins+1;
        }
    }
        
	return s.substr(pos, res);
}

string convert(string s, int nRows) {
    int r, c, col, ins;
    ins = nRows*1.5;
    col = s.length()*2/ins+1;
    char **zigzag;
    string res;

	zigzag = new char* [nRows];
	for ( r=0; r<nRows; r++ ) zigzag[r] = new char [col];

    for ( c=0; c<col; c++ ) {
        for ( r=0; r<nRows; r++ ) {
            if ( c%2==0 && r+c/2*ins<s.length() ) zigzag[r][c] = s[r+c/2*ins];
			else if ( c%2==1 && r%2==1 && nRows+r/2+c/2*ins<s.length() ) zigzag[r][c]=s[nRows+r/2+c/2*ins];
			else zigzag[r][c] = ' ';
        }
    }

    for ( r=0; r<nRows; r++ ) {
        for ( c=0; c<col; c++ )
            if ( zigzag[r][c]!=' ' ) res.push_back(zigzag[r][c]);
    }
        
    return res;
}

int atoi(const char *str) {
    long long value=0;
    int i=-1, symble=0;
        
    while ( str[++i]!='\0' ){
        if ( str[i]==' ' ) continue;
        if ( symble==0 && int(str[i])==45 ) symble = -1;
        if ( symble==0 && int(str[i])==43 ) symble = 1;
    	else if ( symble==0 && int(str[i])>=48 && int(str[i])<=57 ) {
    		value = value*10 + int(str[i]) - 48;
    		symble = 1;
    	}
    		
    	if ( symble==0 ) break;
    
        while ( str[++i]!='\0' && int(str[i])>=48 && int(str[i])<=57 ) value = value*10 + int(str[i]) - 48;
                
        if ( symble!=0 ) break;
    }
        
    value *= symble;
    if ( value>INT_MAX ) return INT_MAX;
    if ( value<INT_MIN ) return INT_MIN;
        
    return value;
}

bool isPalindrome(int x) {
    if ( x<0 ) return false;
    
	int d=pow( 10.0, (int)log10((double)x) );
	while( d>1 ) {
		if ( x%10!=x/d ) return false;
		x = x%d/10;
		d /= 100;
	}
        
    return true;
}

// Regular Expression
bool isMatch(const char *s, const char *p) {
	if (*p == '\0') return *s == '\0';
        
    if (*(p+1) == '*') // next is '*'
    {
        while ((*s == *p || *p == '.') && *s != '\0')
        {
            if (isMatch(s, p+2))
                return true;
            s++;
        }
        return isMatch(s, p+2);
    }
        
    if (*s == '\0') return false;
    return (*s == *p || *p == '.') && isMatch(s+1, p+1);
}

//Sum3, wrong duplicate
vector<vector<int> > threeSum(vector<int> &num) {
    vector< vector<int> > res;
    vector<int> pos, neg, vtmp;
    bool zero = false;
    int tmp, k;
        
    for ( int i=0; i<num.size(); i++ ) {
        if ( num[i]==0 ) zero = true;
        if ( num[i]>0 ) pos.push_back( num[i] );
        if ( num[i]<0 ) neg.push_back( num[i] );
    }
        
    for ( int i=0; i<pos.size(); i++ ) {
        for ( int j=0; j<neg.size(); j++ ) {
            tmp = pos[i]+neg[j];
			vtmp.clear();
            if ( tmp==0 && zero ) {
				vtmp.push_back( neg[j] );
				vtmp.push_back( 0 );
				vtmp.push_back( pos[i] );
				res.push_back( vtmp );
				break;
			}
			else if ( tmp>0 ){
                k=j;
                while ( ++k<neg.size() )
                    if ( tmp+neg[k]==0 ) break;
                if ( k<neg.size() ) {
					vtmp.push_back( min(neg[j], neg[k]) );
					vtmp.push_back( max(neg[j], neg[k]) );
					vtmp.push_back( pos[i] );
                    res.push_back( vtmp );
                }
            }
            else if ( tmp<0 ) {
                k=i;
                while ( ++k<pos.size() )
                    if ( tmp+pos[k]==0 ) break;
                if ( k<pos.size() ) {
					vtmp.push_back( neg[j] );
					vtmp.push_back( min(pos[i], pos[k]) );
					vtmp.push_back( max(pos[i], pos[k]) );
                    res.push_back( vtmp );
                    break;
                }
            }
        }
    }
        
    return res;
}

//Valid Parentheses
char map( const char tmp ){
    if( tmp==')' ) return '(';
	if( tmp=='}' ) return '{';
    if( tmp==']' ) return '[';
        
    return ' ';
}
bool isValid(string s) {
    bool res=true;
    stack<char> brace;
        
    for( int i=0; i<s.length(); ++i ) {
        if( s[i]=='(' || s[i]=='{' || s[i]=='[' ) brace.push( s[i] );
        else if ( s[i]==')' || s[i]=='}' || s[i]==']' ){
			if ( !brace.empty() && map( s[i] )==brace.top() ) brace.pop();
            else return false;
        }
    }
        
    if ( !brace.empty() ) return false;
        
    return res;
}

//
void generateParenthesisRe(int left, int right, string s, vector<string> &res) {
    if (left == 0 && right == 0)
        res.push_back(s);
    if (left > 0)
        generateParenthesisRe(left - 1, right, s + "(", res);
    if (right > left)
        generateParenthesisRe(left, right - 1, s + ")", res);
}

// swap pairs
ListNode *swapPairs(ListNode *head) {
    ListNode list(-1), *tmp, *pos;
        
    list.next = head;
    pos = &list;
        
    if( !(head && head->next)) return head;
        
    while( pos && pos->next && pos->next->next ){
        tmp = pos->next;
        pos->next = pos->next->next;
		tmp->next = pos->next->next;
        pos->next->next = tmp;
            
        pos = tmp;
    }
        
    return list.next;
}

//Reverse Nodes in k-Group
void reverseK(ListNode *prev, int k);
ListNode *reverseKGroup(ListNode *head, int k) {
    ListNode list(-1), *pos;
	int tmp;
        
    list.next = head;
    pos = &list;
        
    while( pos ){
		tmp = k;
        reverseK( pos, tmp);
            
        while( tmp-- && pos ) pos = pos->next;
    }
        
    return list.next;
}
    
void reverseK(ListNode *prev, int k) {
    ListNode *pos, *tmp_1, *tmp_2;
        
    if(k==1) return;
    if( k==2 && prev && prev->next && prev->next->next ) {
        pos = prev->next;
            
        prev->next = pos->next;
        pos->next = prev->next->next;
        prev->next->next = pos;
            
        return;
    }
        
    pos = prev;
	int tmp = k;
    while( --tmp && pos ) pos=pos->next;
    if( pos && pos->next ){
        tmp_1 = prev->next;
            
        prev->next = pos->next;
        pos->next = tmp_1;
            
        tmp_2 = prev->next->next;
        prev->next->next = tmp_1->next;
        tmp_1->next = tmp_2;
            
        reverseK( prev->next, k-2 );
    }
}

//Implement strStr() 
char *strStr(char *haystack, char *needle) {
    char *tmp_h, *tmp_n;
    
	if( haystack==NULL || needle==NULL ) return NULL;
    if( *haystack=='\0' && *needle=='\0' ) return haystack;

    while( *haystack!='\0' ) {
        tmp_n = needle;
        tmp_h = haystack;
        while( *tmp_n!='\0' && *tmp_h!='\0' && *tmp_n==*(tmp_h++) ) tmp_n++;
            
        if( *tmp_n=='\0' ) return haystack;
            
        haystack++;
    }
        
    return NULL;
}


//Divider
int divide(int dividend, int divisor) {
    int result = 0;
    bool flag = ( dividend>0 && divisor<0 || dividend<0 && divisor>0 );
        
    long long dividendll = abs((long long)dividend);
    long long divisorll = abs((long long)divisor);
    while( dividendll>=divisorll ) {
        long long div = divisorll;
        int quat = 1;
            
        while( (div<<1)<=dividendll ) {
            div <<= 1;
            quat <<= 1;
        }
            
        dividendll -= div;
        result += quat;
    }
        
    return flag? -result:result;
}

//Substring with Concatenation of All Words
vector<int> findSubstring(string S, vector<string> &L) {
    vector<int> result;
    if( S.empty() || L.empty() || S.length()<L.size()*L[0].size() ) return result;
    unordered_map<string, int> find, need;
	int len = L[0].size();

	for(auto tmp=L.begin(); tmp!=L.end(); tmp++)
		need[*tmp]++;
        
	for(int i=0; i<=S.length()-L.size()*len; i++) {
		find.clear();
		int j;
		for(j=0; j<L.size(); j++) {
			string sub = S.substr(i+j*len, len);
			auto it = need.find(sub);
			if( it==need.end() ) break;
			if( it->second<=find[sub] ) break;
			find[sub]++;
		}
		if( j==L.size() ) result.push_back(i);
	}
        
    return result;
}

//Longest Valid Parentheses
int longestValidParentheses(string s) {
    int res=0, tmp=0;
    stack<int> match;
        
    for(int i=0; i<s.size(); i++) {
        if( s[i]=='(' ) {
            match.push(tmp);
			tmp=0;
        }
        else if( !match.empty() ) {
			tmp += (1 + match.top());
            match.pop();

			res = max(res, tmp);
        }
        else
			tmp=0;
    }
    
    return res*2;
}

//Search in Rotated Sorted Array
int search(int A[], int n, int target) {
	int i = 0, j = n-1, mid;

	while(i<=j) {
		mid = (i+j)/2;

		if( A[mid]==target ) return mid;
		else if( A[i]<=A[mid] ) {
			if( A[i]<=target && target<A[mid] ) j = mid;
			else i = mid + 1;
		}
		else {
			if( A[i]>target && target>A[mid] ) i = mid + 1;
			else j = mid;
		}
	}

	return -1;
}

//Valid Sudoku
bool isValidSudoku(vector<vector<char> > &board) {
    unordered_map<int, int> check;
        
    if( board.size()!=9 || board[0].size()!=9 ) return false;
        
    //Check for row
    for( int i=0; i<9; i++ ) {
        check.clear();
        for ( int j=0; j<9; j++ ) {
            if( board[i][j]!='.' ) check[board[i][j]-'0']++;
			if( check[board[i][j]-'0']==2 ) return false;
        }
    }
        
    //Check for column
    for( int i=0; i<9; i++ ) {
        check.clear();
        for ( int j=0; j<9; j++ ) {
            if( board[j][i]!='.' ) check[board[j][i]-'0']++;
            if( check[board[j][i]-'0']==2 ) return false;
        }
    }
        
    //Check for 3*3 block
    for( int row=0; row<9; row+=3 )
        for ( int col=0; col<9; col+=3 ) {
			check.clear();
            for (int i=0; i<3; i++) 
				for ( int j=0; j<3; j++ ){
					if( board[row+i][col+j]!='.' ) check[board[row+i][col+j]-'0']++;
					if( check[board[row+i][col+j]-'0']==2 ) return false;
				}
		}
        
        
    return true;
}

//Count and Say
string countAndSay(string now, int n) {
    if( n==0 ) return now;
        
    string next;
    char tmp = now[0];
    char count = '1';
    for( int i=1; i<now.size(); i++ ) {
        if( now[i]==tmp ) count++;
        else {
            next += count; next += tmp;
            tmp = now[i];
            count = '1';
        }
    }
	next += count; next += tmp;
    return countAndSay(next, n-1);
}

string countAndSay(int n) {
    string res;
        
    res = countAndSay("1", n-1);
        
    return res;
}
    
//Trap Water
int trap(int A[], int n) {
    int res=0, tmp;
    int *lmax = new int[n];
    int *rmax = new int[n];
        
    if( n<3 ) return 0;
        
	for(int i=0, tmp=A[0]; i<n; ++i) {
		tmp = max(tmp, A[i]);
		lmax[i] = tmp;
	}

    for(int i=n-1, tmp=A[n-1]; i>=0; --i) {
		tmp = max(tmp, A[i]);
		rmax[i] = tmp;
	}

	for(int i=0; i<n; i++)
		res += min(lmax[i], rmax[i]) - A[i];

	delete [] lmax;
	delete [] rmax;

    return res;
}

//Combination Sum
void combinationSum(const vector<int> &candidates, int pos, int target, vector<int> suspect, vector<vector<int>> &res) {
    if(target==0) {
    	res.push_back( suspect );
    	return;
    }
    	
    for( int i=pos; i<candidates.size() && target>=candidates[i]; i++ ) {
    	suspect.push_back( candidates[i] );
    	combinationSum(candidates, i, target-candidates[i], suspect, res);
    	suspect.pop_back();
    }
}

//Combination Sum2
void combinationSum2(const vector<int> &num, int pos, int target, vector<int> suspect, vector<vector<int>> &res) {
    if(target==0) {
        res.push_back(suspect);
        return;
    }
        
    for(int i=pos; i<num.size() && target>=num[i]; i++) {
        if( i>pos && num[i]==num[i-1] ) continue;
		suspect.push_back(num[i]);
        combinationSum2(num, i+1, target-num[i], suspect, res);
        suspect.pop_back();
    }
        
}

vector<vector<int> > combinationSum(vector<int> &candidates, int target) {
	vector<vector<int>> res;
	if( candidates.size()==0 ) return res;
	sort(candidates.begin(), candidates.end());
	
	vector<int> tmp;
	combinationSum2(candidates, 0, target, tmp, res);

	return res;
}

//First Missing Positive
int firstMissingPositive(int A[], int n) {
    if(n==0) return 1;
        
    int positive=1;
    for(int i=0; i<n; i++)
        if(A[i]>0 && A[i]<n && A[i]!=i+1 && A[i]!=A[A[i]-1]) {
            int tmp = A[i], pos = A[i]-1;
            A[i--] = A[pos];
            A[pos] = tmp;
        }
        
    for(; positive<=n; positive++)
        if(A[positive-1]!=positive) return positive;
        
    return positive;
}

//Multiply String

string multiply(string num1, string num2) {
	int m=num1.size(), n=num2.size();
	string res(m+n,'0');

	int sum, carry;
	for(int i=m-1; i>=0; i--) {
		carry = 0;
		for(int j=n-1; j>=0; j--) {
			sum = carry + (res[i+j+1]-'0') + (num1[i]-'0')*(num2[j]-'0');
			res[i+j+1] = sum%10 + '0';
			carry = sum/10;
		}
		res[i] += carry;
	}

	while(res.size()>1 && res[0]=='0')
		res.erase(res.begin());

	return res;
}

//Paper poker
vector<int> generatePoker(int n) {
	vector<int> res(n, 0);

	queue<int> que;
	for(int i=0; i<n; i++) que.push(i);

	for(int i=0; !que.empty(); i++) {
		int tmp = que.front(); que.pop();
		if( i%2==0 ) res[tmp] = i/2+1;
		else que.push(tmp);
	}

	return res;
}

//Next Permutation
void nextPermutation(vector<int> &num) {
	for(int pos=num.size()-1; pos>0; pos--)
		for(int j=num.size()-1; num[pos-1]<num[pos] && j>=pos; j--)
			if(num[j]>num[pos-1]) {
				swap(num[j], num[pos-1]);
				sort(num.begin()+pos, num.end());
				return;
			}
    
    sort(num.begin(), num.end());
}

//Anagrams
vector<string> anagrams(vector<string> &strs) {
    vector<string> res;
    unordered_map<string, vector<int>> record;

	for(int i=0; i<strs.size(); i++) {
		string s = strs[i];
		sort(s.begin(), s.end());
		record[s].push_back(i);
	}

	for(auto tmp=record.begin(); tmp!=record.end(); tmp++) {
		vector<int> group = tmp->second;
		for(int i=0; group.size()>1 && i<group.size(); i++)
			res.push_back(strs[group[i]]);
	}

    return res;
}


int main()
{
	char stop;
	Solution solution;
	
	/*int g1[] = { -361, -425, -367, 381, -264, 473, 411, -218, -376, -74, -83, 329, -367, 313, -397, 402, -245, -437, -177, -453, 324, 142, -319, 160, 16, 488, -297, 120, -156, 489, 91, 325, 115, 180, 50, -193, 230, 424, 198, -75, 333, -408, 425, -103, -460, -188, -43, 268, 302, -173, 186 };
	int g2[] = { -314, 199, 363, 116, 325, -331, 79, 189, -317, -156, 297, -66, 282, 272, -292, 35, -84, 10, -217, -130, -232, -283, -61, 168, 418, 146, 247, 343, -234, 95, 463, 5, -68, -367, -371, -257, -355, 244, 114, -134, 311, -443, -448, -461, 402, -472, 1, -279, -315, -30, 204 };
	int g3[] = { -420, 166, -299, 43, 317, -415, 139, -482, 29, 375, 188, 269, -389, 416, 73, -470, 110, 88, -168, -40, -310, 273, -83, -348, -175, -332, -376, 383, 76, 428, 53, 165, -136, -263, 443, 493, 328, 292, 366, 79, -12, -352, -368, 450, 366, -40, 92, 484, 7, -479, -292 };
	int g4[] = { -255, 142, 240, 249, -492, -153, -108, -446, 37, 367, 0, -363, -35, 415, 299, -180, 240, 269, -487, -120, 236, -188, 344, -359, 98, -250, -463, 18, -274, -415, 28, 72, -190, 326, 376, 184, -256, 200, 462, -284, -328, 85, 417, 211, 31, -232, -141, 352, 485, 379, -303 };
	int g5[] = { 197, -250, 9, 131, 268, -320, -249, 426, 477, -466, -125, -25, 358, -154, 436, 93, -381, -150, 269, -466, 275, -279, 397, 191, 45, -339, 418, 135, -475, 20, -371, 77, -240, 370, -449, 98, -433, 407, -446, 27, -322, -250, -394, -408, 150, 204, -219, 21, 315, 370, 302 };
	
	vector<vector<int>> grid;
	grid.push_back(vector<int>(g1, g1 + sizeof(g1) / sizeof(int)));
	grid.push_back(vector<int>(g2, g2 + sizeof(g2) / sizeof(int)));
	grid.push_back(vector<int>(g3, g3 + sizeof(g3) / sizeof(int)));
	grid.push_back(vector<int>(g4, g4 + sizeof(g4) / sizeof(int)));
	grid.push_back(vector<int>(g5, g5 + sizeof(g5) / sizeof(int)));
	*/

	//string theSet[] = { "" };
	//vector<string> s(theSet, theSet + sizeof(theSet)/sizeof(string));
	//vector<string> res = solution.fullJustify(s, 30);
	/*for (int i = 0; i < res.size(); i++) {
		for (int j = 0; j < res[i].size(); j++)
			cout << res[i][j] << '\t';
		cout << endl;
	}*/
	cout << solution.uniquePaths(3, 7) << endl;
	
	cin>>stop;
	return 0;
}