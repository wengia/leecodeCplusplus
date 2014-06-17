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

	string start = "nanny", end = "aloud";
	string theSet[] = { "ricky", "grind", "cubic", "panic", "lover", "farce", "gofer", "sales", "flint", "omens", "lipid", "briny", "cloth", "anted", "slime", "oaten", "harsh", "touts", "stoop", "cabal", "lazed", "elton", "skunk", "nicer", "pesky", "kusch", "bused", "kinda", "tunis", "enjoy", "aches", "prowl", "babar", "rooms", "burst", "slush", "pines", "urine", "pinky", "bayed", "mania", "light", "flare", "wares", "women", "verne", "moron", "shine", "bluer", "zeros", "bleak", "brief", "tamra", "vasts", "jamie", "lairs", "penal", "worst", "yowls", "pills", "taros", "addle", "alyce", "creep", "saber", "floyd", "cures", "soggy", "vexed", "vilma", "cabby", "verde", "euler", "cling", "wanna", "jenny", "donor", "stole", "sakha", "blake", "sanes", "riffs", "forge", "horus", "sered", "piked", "prosy", "wases", "glove", "onset", "spake", "benin", "talks", "sites", "biers", "wendy", "dante", "allan", "haven", "nears", "shaka", "sloth", "perky", "spear", "spend", "clint", "dears", "sadly", "units", "vista", "hinds", "marat", "natal", "least", "bough", "pales", "boole", "ditch", "greys", "slunk", "bitch", "belts", "sense", "skits", "monty", "yawns", "music", "hails", "alien", "gibes", "lille", "spacy", "argot", "wasps", "drubs", "poops", "bella", "clone", "beast", "emend", "iring", "start", "darla", "bells", "cults", "dhaka", "sniff", "seers", "bantu", "pages", "fever", "tacky", "hoses", "strop", "climb", "pairs", "later", "grant", "raven", "stael", "drips", "lucid", "awing", "dines", "balms", "della", "galen", "toned", "snips", "shady", "chili", "fears", "nurse", "joint", "plump", "micky", "lions", "jamal", "queer", "ruins", "frats", "spoof", "semen", "pulps", "oldie", "coors", "rhone", "papal", "seals", "spans", "scaly", "sieve", "klaus", "drums", "tided", "needs", "rider", "lures", "treks", "hares", "liner", "hokey", "boots", "primp", "laval", "limes", "putts", "fonda", "damon", "pikes", "hobbs", "specs", "greet", "ketch", "braid", "purer", "tsars", "berne", "tarts", "clean", "grate", "trips", "chefs", "timex", "vicky", "pares", "price", "every", "beret", "vices", "jodie", "fanny", "mails", "built", "bossy", "farms", "pubic", "gongs", "magma", "quads", "shell", "jocks", "woods", "waded", "parka", "jells", "worse", "diner", "risks", "bliss", "bryan", "terse", "crier", "incur", "murky", "gamed", "edges", "keens", "bread", "raced", "vetch", "glint", "zions", "porno", "sizes", "mends", "ached", "allie", "bands", "plank", "forth", "fuels", "rhyme", "wimpy", "peels", "foggy", "wings", "frill", "edgar", "slave", "lotus", "point", "hints", "germs", "clung", "limed", "loafs", "realm", "myron", "loopy", "plush", "volts", "bimbo", "smash", "windy", "sours", "choke", "karin", "boast", "whirr", "tiber", "dimes", "basel", "cutes", "pinto", "troll", "thumb", "decor", "craft", "tared", "split", "josue", "tramp", "screw", "label", "lenny", "apses", "slept", "sikhs", "child", "bouts", "cites", "swipe", "lurks", "seeds", "fists", "hoard", "steed", "reams", "spoil", "diego", "peale", "bevel", "flags", "mazes", "quart", "snipe", "latch", "lards", "acted", "falls", "busby", "holed", "mummy", "wrong", "wipes", "carlo", "leers", "wails", "night", "pasty", "eater", "flunk", "vedas", "curse", "tyros", "mirth", "jacky", "butte", "wired", "fixes", "tares", "vague", "roved", "stove", "swoon", "scour", "coked", "marge", "cants", "comic", "corns", "zilch", "typos", "lives", "truer", "comma", "gaily", "teals", "witty", "hyper", "croat", "sways", "tills", "hones", "dowel", "llano", "clefs", "fores", "cinch", "brock", "vichy", "bleed", "nuder", "hoyle", "slams", "macro", "arabs", "tauts", "eager", "croak", "scoop", "crime", "lurch", "weals", "fates", "clipt", "teens", "bulls", "domed", "ghana", "culls", "frame", "hanky", "jared", "swain", "truss", "drank", "lobby", "lumps", "pansy", "whews", "saris", "trite", "weeps", "dozes", "jeans", "flood", "chimu", "foxes", "gelds", "sects", "scoff", "poses", "mares", "famed", "peers", "hells", "laked", "zests", "wring", "steal", "snoot", "yodel", "scamp", "ellis", "bandy", "marry", "jives", "vises", "blurb", "relay", "patch", "haley", "cubit", "heine", "place", "touch", "grain", "gerry", "badly", "hooke", "fuchs", "savor", "apron", "judge", "loren", "britt", "smith", "tammy", "altar", "duels", "huber", "baton", "dived", "apace", "sedan", "basts", "clark", "mired", "perch", "hulks", "jolly", "welts", "quack", "spore", "alums", "shave", "singe", "lanny", "dread", "profs", "skeet", "flout", "darin", "newed", "steer", "taine", "salvo", "mites", "rules", "crash", "thorn", "olive", "saves", "yawed", "pique", "salon", "ovens", "dusty", "janie", "elise", "carve", "winds", "abash", "cheep", "strap", "fared", "discs", "poxed", "hoots", "catch", "combo", "maize", "repay", "mario", "snuff", "delve", "cored", "bards", "sudan", "shuns", "yukon", "jowls", "wayne", "torus", "gales", "creek", "prove", "needy", "wisps", "terri", "ranks", "books", "dicky", "tapes", "aping", "padre", "roads", "nines", "seats", "flats", "rains", "moira", "basic", "loves", "pulls", "tough", "gills", "codes", "chest", "teeny", "jolts", "woody", "flame", "asked", "dulls", "hotly", "glare", "mucky", "spite", "flake", "vines", "lindy", "butts", "froth", "beeps", "sills", "bunny", "flied", "shaun", "mawed", "velds", "voled", "doily", "patel", "snake", "thigh", "adler", "calks", "desks", "janus", "spunk", "baled", "match", "strip", "hosed", "nippy", "wrest", "whams", "calfs", "sleet", "wives", "boars", "chain", "table", "duked", "riped", "edens", "galas", "huffs", "biddy", "claps", "aleut", "yucks", "bangs", "quids", "glenn", "evert", "drunk", "lusts", "senna", "slate", "manet", "roted", "sleep", "loxes", "fluky", "fence", "clamp", "doted", "broad", "sager", "spark", "belch", "mandy", "deana", "beyer", "hoist", "leafy", "levee", "libel", "tonic", "aloes", "steam", "skews", "tides", "stall", "rifts", "saxon", "mavis", "asama", "might", "dotes", "tangs", "wroth", "kited", "salad", "liens", "clink", "glows", "balky", "taffy", "sided", "sworn", "oasis", "tenth", "blurt", "tower", "often", "walsh", "sonny", "andes", "slump", "scans", "boded", "chive", "finer", "ponce", "prune", "sloes", "dined", "chums", "dingo", "harte", "ahead", "event", "freer", "heart", "fetch", "sated", "soapy", "skins", "royal", "cuter", "loire", "minot", "aisle", "horny", "slued", "panel", "eight", "snoop", "pries", "clive", "pored", "wrist", "piped", "daren", "cells", "parks", "slugs", "cubed", "highs", "booze", "weary", "stain", "hoped", "finny", "weeds", "fetid", "racer", "tasks", "right", "saint", "shahs", "basis", "refer", "chart", "seize", "lulls", "slant", "belay", "clots", "jinny", "tours", "modes", "gloat", "dunks", "flute", "conch", "marts", "aglow", "gayer", "lazes", "dicks", "chime", "bears", "sharp", "hatch", "forms", "terry", "gouda", "thins", "janet", "tonya", "axons", "sewed", "danny", "rowdy", "dolts", "hurry", "opine", "fifty", "noisy", "spiky", "humid", "verna", "poles", "jayne", "pecos", "hooky", "haney", "shams", "snots", "sally", "ruder", "tempe", "plunk", "shaft", "scows", "essie", "dated", "fleet", "spate", "bunin", "hikes", "sodas", "filly", "thyme", "fiefs", "perks", "chary", "kiths", "lidia", "lefty", "wolff", "withe", "three", "crawl", "wotan", "brown", "japed", "tolls", "taken", "threw", "crave", "clash", "layer", "tends", "notes", "fudge", "musky", "bawdy", "aline", "matts", "shirr", "balks", "stash", "wicks", "crepe", "foods", "fares", "rotes", "party", "petty", "press", "dolly", "mangy", "leeks", "silly", "leant", "nooks", "chapt", "loose", "caged", "wages", "grist", "alert", "sheri", "moody", "tamps", "hefts", "souls", "rubes", "rolex", "skulk", "veeps", "nonce", "state", "level", "whirl", "bight", "grits", "reset", "faked", "spiny", "mixes", "hunks", "major", "missy", "arius", "damns", "fitly", "caped", "mucus", "trace", "surat", "lloyd", "furry", "colin", "texts", "livia", "reply", "twill", "ships", "peons", "shear", "norms", "jumbo", "bring", "masks", "zippy", "brine", "dorks", "roded", "sinks", "river", "wolfs", "strew", "myths", "pulpy", "prank", "veins", "flues", "minus", "phone", "banns", "spell", "burro", "brags", "boyle", "lambs", "sides", "knees", "clews", "aired", "skirt", "heavy", "dimer", "bombs", "scums", "hayes", "chaps", "snugs", "dusky", "loxed", "ellen", "while", "swank", "track", "minim", "wiled", "hazed", "roofs", "cantu", "sorry", "roach", "loser", "brass", "stint", "jerks", "dirks", "emory", "campy", "poise", "sexed", "gamer", "catty", "comte", "bilbo", "fasts", "ledge", "drier", "idles", "doors", "waged", "rizal", "pured", "weirs", "crisp", "tasty", "sored", "palmy", "parts", "ethel", "unify", "crows", "crest", "udder", "delis", "punks", "dowse", "totes", "emile", "coded", "shops", "poppa", "pours", "gushy", "tiffs", "shads", "birds", "coils", "areas", "boons", "hulls", "alter", "lobes", "pleat", "depth", "fires", "pones", "serra", "sweat", "kline", "malay", "ruled", "calve", "tired", "drabs", "tubed", "wryer", "slung", "union", "sonya", "aided", "hewed", "dicey", "grids", "nixed", "whits", "mills", "buffs", "yucky", "drops", "ready", "yuppy", "tweet", "napes", "cadre", "teach", "rasps", "dowdy", "hoary", "canto", "posed", "dumbo", "kooks", "reese", "snaky", "binge", "byron", "phony", "safer", "friar", "novel", "scale", "huron", "adorn", "carla", "fauna", "myers", "hobby", "purse", "flesh", "smock", "along", "boils", "pails", "times", "panza", "lodge", "clubs", "colby", "great", "thing", "peaks", "diana", "vance", "whets", "bergs", "sling", "spade", "soaks", "beach", "traps", "aspen", "romps", "boxed", "fakir", "weave", "nerds", "swazi", "dotty", "curls", "diver", "jonas", "waite", "verbs", "yeast", "lapel", "barth", "soars", "hooks", "taxed", "slews", "gouge", "slags", "chang", "chafe", "saved", "josie", "syncs", "fonds", "anion", "actor", "seems", "pyrex", "isiah", "glued", "groin", "goren", "waxes", "tonia", "whine", "scads", "knelt", "teaks", "satan", "tromp", "spats", "merry", "wordy", "stake", "gland", "canal", "donna", "lends", "filed", "sacks", "shied", "moors", "paths", "older", "pooch", "balsa", "riced", "facet", "decaf", "attic", "elder", "akron", "chomp", "chump", "picky", "money", "sheer", "bolls", "crabs", "dorms", "water", "veers", "tease", "dummy", "dumbs", "lethe", "halls", "rifer", "demon", "fucks", "whips", "plops", "fuses", "focal", "taces", "snout", "edict", "flush", "burps", "dawes", "lorry", "spews", "sprat", "click", "deann", "sited", "aunts", "quips", "godly", "pupil", "nanny", "funks", "shoon", "aimed", "stacy", "helms", "mints", "banks", "pinch", "local", "twine", "pacts", "deers", "halos", "slink", "preys", "potty", "ruffs", "pusan", "suits", "finks", "slash", "prods", "dense", "edsel", "heeds", "palls", "slats", "snits", "mower", "rares", "ailed", "rouge", "ellie", "gated", "lyons", "duded", "links", "oaths", "letha", "kicks", "firms", "gravy", "month", "kongo", "mused", "ducal", "toted", "vocal", "disks", "spied", "studs", "macao", "erick", "coupe", "starr", "reaps", "decoy", "rayon", "nicks", "breed", "cosby", "haunt", "typed", "plain", "trays", "muled", "saith", "drano", "cower", "snows", "buses", "jewry", "argus", "doers", "flays", "swish", "resin", "boobs", "sicks", "spies", "bails", "wowed", "mabel", "check", "vapid", "bacon", "wilda", "ollie", "loony", "irked", "fraud", "doles", "facts", "lists", "gazed", "furls", "sunks", "stows", "wilde", "brick", "bowed", "guise", "suing", "gates", "niter", "heros", "hyped", "clomp", "never", "lolls", "rangy", "paddy", "chant", "casts", "terns", "tunas", "poker", "scary", "maims", "saran", "devon", "tripe", "lingo", "paler", "coped", "bride", "voted", "dodge", "gross", "curds", "sames", "those", "tithe", "steep", "flaks", "close", "swops", "stare", "notch", "prays", "roles", "crush", "feuds", "nudge", "baned", "brake", "plans", "weepy", "dazed", "jenna", "weiss", "tomes", "stews", "whist", "gibed", "death", "clank", "cover", "peeks", "quick", "abler", "daddy", "calls", "scald", "lilia", "flask", "cheer", "grabs", "megan", "canes", "jules", "blots", "mossy", "begun", "freak", "caved", "hello", "hades", "theed", "wards", "darcy", "malta", "peter", "whorl", "break", "downs", "odder", "hoofs", "kiddo", "macho", "fords", "liked", "flees", "swing", "elect", "hoods", "pluck", "brook", "astir", "bland", "sward", "modal", "flown", "ahmad", "waled", "craps", "cools", "roods", "hided", "plath", "kings", "grips", "gives", "gnats", "tabby", "gauls", "think", "bully", "fogey", "sawed", "lints", "pushy", "banes", "drake", "trail", "moral", "daley", "balds", "chugs", "geeky", "darts", "soddy", "haves", "opens", "rends", "buggy", "moles", "freud", "gored", "shock", "angus", "puree", "raves", "johns", "armed", "packs", "minis", "reich", "slots", "totem", "clown", "popes", "brute", "hedge", "latin", "stoke", "blend", "pease", "rubik", "greer", "hindi", "betsy", "flows", "funky", "kelli", "humps", "chewy", "welds", "scowl", "yells", "cough", "sasha", "sheaf", "jokes", "coast", "words", "irate", "hales", "camry", "spits", "burma", "rhine", "bends", "spill", "stubs", "power", "voles", "learn", "knoll", "style", "twila", "drove", "dacca", "sheen", "papas", "shale", "jones", "duped", "tunny", "mouse", "floss", "corks", "skims", "swaps", "inned", "boxer", "synch", "skies", "strep", "bucks", "belau", "lower", "flaky", "quill", "aural", "rufus", "floes", "pokes", "sends", "sates", "dally", "boyer", "hurts", "foyer", "gowns", "torch", "luria", "fangs", "moats", "heinz", "bolts", "filet", "firth", "begot", "argue", "youth", "chimp", "frogs", "kraft", "smite", "loges", "loons", "spine", "domes", "pokey", "timur", "noddy", "doggy", "wades", "lanes", "hence", "louts", "turks", "lurid", "goths", "moist", "bated", "giles", "stood", "winos", "shins", "potts", "brant", "vised", "alice", "rosie", "dents", "babes", "softy", "decay", "meats", "tanya", "rusks", "pasts", "karat", "nuked", "gorge", "kinks", "skull", "noyce", "aimee", "watch", "cleat", "stuck", "china", "testy", "doses", "safes", "stage", "bayes", "twins", "limps", "denis", "chars", "flaps", "paces", "abase", "grays", "deans", "maria", "asset", "smuts", "serbs", "whigs", "vases", "robyn", "girls", "pents", "alike", "nodal", "molly", "swigs", "swill", "slums", "rajah", "bleep", "beget", "thanh", "finns", "clock", "wafts", "wafer", "spicy", "sorer", "reach", "beats", "baker", "crown", "drugs", "daisy", "mocks", "scots", "fests", "newer", "agate", "drift", "marta", "chino", "flirt", "homed", "bribe", "scram", "bulks", "servo", "vesta", "divas", "preps", "naval", "tally", "shove", "ragas", "blown", "droll", "tryst", "lucky", "leech", "lines", "sires", "pyxed", "taper", "trump", "payee", "midge", "paris", "bored", "loads", "shuts", "lived", "swath", "snare", "boned", "scars", "aeons", "grime", "writs", "paige", "rungs", "blent", "signs", "davis", "dials", "daubs", "rainy", "fawns", "wrier", "golds", "wrath", "ducks", "allow", "hosea", "spike", "meals", "haber", "muses", "timed", "broom", "burks", "louis", "gangs", "pools", "vales", "altai", "elope", "plied", "slain", "chasm", "entry", "slide", "bawls", "title", "sings", "grief", "viola", "doyle", "peach", "davit", "bench", "devil", "latex", "miles", "pasha", "tokes", "coves", "wheel", "tried", "verdi", "wanda", "sivan", "prior", "fryer", "plots", "kicky", "porch", "shill", "coats", "borne", "brink", "pawed", "erwin", "tense", "stirs", "wends", "waxen", "carts", "smear", "rival", "scare", "phase", "bragg", "crane", "hocks", "conan", "bests", "dares", "molls", "roots", "dunes", "slips", "waked", "fours", "bolds", "slosh", "yemen", "poole", "solid", "ports", "fades", "legal", "cedes", "green", "curie", "seedy", "riper", "poled", "glade", "hosts", "tools", "razes", "tarry", "muddy", "shims", "sword", "thine", "lasts", "bloat", "soled", "tardy", "foots", "skiff", "volta", "murks", "croci", "gooks", "gamey", "pyxes", "poems", "kayla", "larva", "slaps", "abuse", "pings", "plows", "geese", "minks", "derby", "super", "inked", "manic", "leaks", "flops", "lajos", "fuzes", "swabs", "twigs", "gummy", "pyres", "shrew", "islet", "doled", "wooly", "lefts", "hunts", "toast", "faith", "macaw", "sonia", "leafs", "colas", "conks", "altos", "wiped", "scene", "boors", "patsy", "meany", "chung", "wakes", "clear", "ropes", "tahoe", "zones", "crate", "tombs", "nouns", "garth", "puked", "chats", "hanks", "baked", "binds", "fully", "soaps", "newel", "yarns", "puers", "carps", "spelt", "lully", "towed", "scabs", "prime", "blest", "patty", "silky", "abner", "temps", "lakes", "tests", "alias", "mines", "chips", "funds", "caret", "splat", "perry", "turds", "junks", "cramp", "saned", "peary", "snarl", "fired", "stung", "nancy", "bulge", "styli", "seams", "hived", "feast", "triad", "jaded", "elvin", "canny", "birth", "routs", "rimed", "pusey", "laces", "taste", "basie", "malls", "shout", "prier", "prone", "finis", "claus", "loops", "heron", "frump", "spare", "menus", "ariel", "crams", "bloom", "foxed", "moons", "mince", "mixed", "piers", "deres", "tempt", "dryer", "atone", "heats", "dario", "hawed", "swims", "sheet", "tasha", "dings", "clare", "aging", "daffy", "wried", "foals", "lunar", "havel", "irony", "ronny", "naves", "selma", "gurus", "crust", "percy", "murat", "mauro", "cowed", "clang", "biker", "harms", "barry", "thump", "crude", "ulnae", "thong", "pager", "oases", "mered", "locke", "merle", "soave", "petal", "poser", "store", "winch", "wedge", "inlet", "nerdy", "utter", "filth", "spray", "drape", "pukes", "ewers", "kinds", "dates", "meier", "tammi", "spoor", "curly", "chill", "loped", "gooey", "boles", "genet", "boost", "beets", "heath", "feeds", "growl", "livid", "midst", "rinds", "fresh", "waxed", "yearn", "keeps", "rimes", "naked", "flick", "plies", "deeps", "dirty", "hefty", "messy", "hairy", "walks", "leper", "sykes", "nerve", "rover", "jived", "brisk", "lenin", "viper", "chuck", "sinus", "luger", "ricks", "hying", "rusty", "kathy", "herds", "wider", "getty", "roman", "sandy", "pends", "fezes", "trios", "bites", "pants", "bless", "diced", "earth", "shack", "hinge", "melds", "jonah", "chose", "liver", "salts", "ratty", "ashed", "wacky", "yokes", "wanly", "bruce", "vowel", "black", "grail", "lungs", "arise", "gluts", "gluey", "navel", "coyer", "ramps", "miter", "aldan", "booth", "musty", "rills", "darns", "tined", "straw", "kerri", "hared", "lucks", "metes", "penny", "radon", "palms", "deeds", "earls", "shard", "pried", "tampa", "blank", "gybes", "vicki", "drool", "groom", "curer", "cubes", "riggs", "lanky", "tuber", "caves", "acing", "golly", "hodge", "beard", "ginny", "jibed", "fumes", "astor", "quito", "cargo", "randi", "gawky", "zings", "blind", "dhoti", "sneak", "fatah", "fixer", "lapps", "cline", "grimm", "fakes", "maine", "erika", "dealt", "mitch", "olden", "joist", "gents", "likes", "shelf", "silts", "goats", "leads", "marin", "spire", "louie", "evans", "amuse", "belly", "nails", "snead", "model", "whats", "shari", "quote", "tacks", "nutty", "lames", "caste", "hexes", "cooks", "miner", "shawn", "anise", "drama", "trike", "prate", "ayers", "loans", "botch", "vests", "cilia", "ridge", "thugs", "outed", "jails", "moped", "plead", "tunes", "nosed", "wills", "lager", "lacks", "cried", "wince", "berle", "flaws", "boise", "tibet", "bided", "shred", "cocky", "brice", "delta", "congo", "holly", "hicks", "wraps", "cocks", "aisha", "heard", "cured", "sades", "horsy", "umped", "trice", "dorky", "curve", "ferry", "haler", "ninth", "pasta", "jason", "honer", "kevin", "males", "fowls", "awake", "pores", "meter", "skate", "drink", "pussy", "soups", "bases", "noyes", "torts", "bogus", "still", "soupy", "dance", "worry", "eldon", "stern", "menes", "dolls", "dumpy", "gaunt", "grove", "coops", "mules", "berry", "sower", "roams", "brawl", "greed", "stags", "blurs", "swift", "treed", "taney", "shame", "easel", "moves", "leger", "ville", "order", "spock", "nifty", "brian", "elias", "idler", "serve", "ashen", "bizet", "gilts", "spook", "eaten", "pumas", "cotes", "broke", "toxin", "groan", "laths", "joins", "spots", "hated", "tokay", "elite", "rawer", "fiats", "cards", "sassy", "milks", "roost", "glean", "lutes", "chins", "drown", "marks", "pined", "grace", "fifth", "lodes", "rusts", "terms", "maxes", "savvy", "choir", "savoy", "spoon", "halve", "chord", "hulas", "sarah", "celia", "deems", "ninny", "wines", "boggy", "birch", "raved", "wales", "beams", "vibes", "riots", "warty", "nigel", "askew", "faxes", "sedge", "sheol", "pucks", "cynic", "relax", "boers", "whims", "bents", "candy", "luann", "slogs", "bonny", "barns", "iambs", "fused", "duffy", "guilt", "bruin", "pawls", "penis", "poppy", "owing", "tribe", "tuner", "moray", "timid", "ceded", "geeks", "kites", "curio", "puffy", "perot", "caddy", "peeve", "cause", "dills", "gavel", "manse", "joker", "lynch", "crank", "golda", "waits", "wises", "hasty", "paves", "grown", "reedy", "crypt", "tonne", "jerky", "axing", "swept", "posse", "rings", "staff", "tansy", "pared", "glaze", "grebe", "gonna", "shark", "jumps", "vials", "unset", "hires", "tying", "lured", "motes", "linen", "locks", "mamas", "nasty", "mamie", "clout", "nader", "velma", "abate", "tight", "dales", "serer", "rives", "bales", "loamy", "warps", "plato", "hooch", "togae", "damps", "ofter", "plumb", "fifes", "filmy", "wiper", "chess", "lousy", "sails", "brahe", "ounce", "flits", "hindu", "manly", "beaux", "mimed", "liken", "forts", "jambs", "peeps", "lelia", "brews", "handy", "lusty", "brads", "marne", "pesos", "earle", "arson", "scout", "showy", "chile", "sumps", "hiked", "crook", "herbs", "silks", "alamo", "mores", "dunce", "blaze", "stank", "haste", "howls", "trots", "creon", "lisle", "pause", "hates", "mulch", "mined", "moder", "devin", "types", "cindy", "beech", "tuned", "mowed", "pitts", "chaos", "colds", "bidet", "tines", "sighs", "slimy", "brain", "belle", "leery", "morse", "ruben", "prows", "frown", "disco", "regal", "oaken", "sheds", "hives", "corny", "baser", "fated", "throe", "revel", "bores", "waved", "shits", "elvia", "ferns", "maids", "color", "coifs", "cohan", "draft", "hmong", "alton", "stine", "cluck", "nodes", "emily", "brave", "blair", "blued", "dress", "bunts", "holst", "clogs", "rally", "knack", "demos", "brady", "blues", "flash", "goofy", "blocs", "diane", "colic", "smile", "yules", "foamy", "splay", "bilge", "faker", "foils", "condo", "knell", "crack", "gallo", "purls", "auras", "cakes", "doves", "joust", "aides", "lades", "muggy", "tanks", "middy", "tarps", "slack", "capet", "frays", "donny", "venal", "yeats", "misty", "denim", "glass", "nudes", "seeps", "gibbs", "blows", "bobbi", "shane", "yards", "pimps", "clued", "quiet", "witch", "boxes", "prawn", "kerry", "torah", "kinko", "dingy", "emote", "honor", "jelly", "grins", "trope", "vined", "bagel", "arden", "rapid", "paged", "loved", "agape", "mural", "budge", "ticks", "suers", "wendi", "slice", "salve", "robin", "bleat", "batik", "myles", "teddy", "flatt", "puppy", "gelid", "largo", "attar", "polls", "glide", "serum", "fundy", "sucks", "shalt", "sewer", "wreak", "dames", "fonts", "toxic", "hines", "wormy", "grass", "louse", "bowls", "crass", "benny", "moire", "margo", "golfs", "smart", "roxie", "wight", "reign", "dairy", "clops", "paled", "oddly", "sappy", "flair", "shown", "bulgy", "benet", "larch", "curry", "gulfs", "fends", "lunch", "dukes", "doris", "spoke", "coins", "manna", "conga", "jinns", "eases", "dunno", "tisha", "swore", "rhino", "calms", "irvin", "clans", "gully", "liege", "mains", "besot", "serge", "being", "welch", "wombs", "draco", "lynda", "forty", "mumps", "bloch", "ogden", "knits", "fussy", "alder", "danes", "loyal", "valet", "wooer", "quire", "liefs", "shana", "toyed", "forks", "gages", "slims", "cloys", "yates", "rails", "sheep", "nacho", "divan", "honks", "stone", "snack", "added", "basal", "hasps", "focus", "alone", "laxes", "arose", "lamed", "wrapt", "frail", "clams", "plait", "hover", "tacos", "mooch", "fault", "teeth", "marva", "mucks", "tread", "waves", "purim", "boron", "horde", "smack", "bongo", "monte", "swirl", "deals", "mikes", "scold", "muter", "sties", "lawns", "fluke", "jilts", "meuse", "fives", "sulky", "molds", "snore", "timmy", "ditty", "gasps", "kills", "carey", "jawed", "byers", "tommy", "homer", "hexed", "dumas", "given", "mewls", "smelt", "weird", "speck", "merck", "keats", "draws", "trent", "agave", "wells", "chews", "blabs", "roves", "grieg", "evens", "alive", "mulls", "cared", "garbo", "fined", "happy", "trued", "rodes", "thurs", "cadet", "alvin", "busch", "moths", "guild", "staci", "lever", "widen", "props", "hussy", "lamer", "riley", "bauer", "chirp", "rants", "poxes", "shyer", "pelts", "funny", "slits", "tinge", "ramos", "shift", "caper", "credo", "renal", "veils", "covey", "elmer", "mated", "tykes", "wooed", "briar", "gears", "foley", "shoes", "decry", "hypes", "dells", "wilds", "runts", "wilts", "white", "easts", "comer", "sammy", "lochs", "favor", "lance", "dawns", "bushy", "muted", "elsie", "creel", "pocks", "tenet", "cagey", "rides", "socks", "ogled", "soils", "sofas", "janna", "exile", "barks", "frank", "takes", "zooms", "hakes", "sagan", "scull", "heaps", "augur", "pouch", "blare", "bulbs", "wryly", "homey", "tubas", "limbo", "hardy", "hoagy", "minds", "bared", "gabby", "bilks", "float", "limns", "clasp", "laura", "range", "brush", "tummy", "kilts", "cooed", "worms", "leary", "feats", "robes", "suite", "veals", "bosch", "moans", "dozen", "rarer", "slyer", "cabin", "craze", "sweet", "talon", "treat", "yanks", "react", "creed", "eliza", "sluts", "cruet", "hafts", "noise", "seder", "flies", "weeks", "venus", "backs", "eider", "uriel", "vouch", "robed", "hacks", "perth", "shiny", "stilt", "torte", "throb", "merer", "twits", "reeds", "shawl", "clara", "slurs", "mixer", "newts", "fried", "woolf", "swoop", "kaaba", "oozed", "mayer", "caned", "laius", "lunge", "chits", "kenny", "lifts", "mafia", "sowed", "piled", "stein", "whack", "colts", "warms", "cleft", "girds", "seeks", "poets", "angel", "trade", "parsi", "tiers", "rojas", "vexes", "bryce", "moots", "grunt", "drain", "lumpy", "stabs", "poohs", "leapt", "polly", "cuffs", "giddy", "towns", "dacha", "quoth", "provo", "dilly", "carly", "mewed", "tzars", "crock", "toked", "speak", "mayas", "pssts", "ocher", "motel", "vogue", "camps", "tharp", "taunt", "drone", "taint", "badge", "scott", "scats", "bakes", "antes", "gruel", "snort", "capes", "plate", "folly", "adobe", "yours", "papaw", "hench", "moods", "clunk", "chevy", "tomas", "narcs", "vonda", "wiles", "prigs", "chock", "laser", "viced", "stiff", "rouse", "helps", "knead", "gazer", "blade", "tumid", "avail", "anger", "egged", "guide", "goads", "rabin", "toddy", "gulps", "flank", "brats", "pedal", "junky", "marco", "tinny", "tires", "flier", "satin", "darth", "paley", "gumbo", "rared", "muffs", "rower", "prude", "frees", "quays", "homes", "munch", "beefs", "leash", "aston", "colon", "finch", "bogey", "leaps", "tempo", "posts", "lined", "gapes", "locus", "maori", "nixes", "liven", "songs", "opted", "babel", "wader", "barer", "farts", "lisps", "koran", "lathe", "trill", "smirk", "mamma", "viler", "scurf", "ravel", "brigs", "cooky", "sachs", "fulls", "goals", "turfs", "norse", "hauls", "cores", "fairy", "pluto", "kneed", "cheek", "pangs", "risen", "czars", "milne", "cribs", "genes", "wefts", "vents", "sages", "seres", "owens", "wiley", "flume", "haded", "auger", "tatty", "onion", "cater", "wolfe", "magic", "bodes", "gulls", "gazes", "dandy", "snags", "rowed", "quell", "spurn", "shore", "veldt", "turns", "slavs", "coach", "stalk", "snuck", "piles", "orate", "joyed", "daily", "crone", "wager", "solos", "earns", "stark", "lauds", "kasey", "villa", "gnaws", "scent", "wears", "fains", "laced", "tamer", "pipes", "plant", "lorie", "rivet", "tamed", "cozen", "theme", "lifer", "sunny", "shags", "flack", "gassy", "eased", "jeeps", "shire", "fargo", "timer", "brash", "behan", "basin", "volga", "krone", "swiss", "docks", "booed", "ebert", "gusty", "delay", "oared", "grady", "buick", "curbs", "crete", "lucas", "strum", "besom", "gorse", "troth", "donne", "chink", "faced", "ahmed", "texas", "longs", "aloud", "bethe", "cacao", "hilda", "eagle", "karyn", "harks", "adder", "verse", "drays", "cello", "taped", "snide", "taxis", "kinky", "penes", "wicca", "sonja", "aways", "dyers", "bolas", "elfin", "slope", "lamps", "hutch", "lobed", "baaed", "masts", "ashes", "ionic", "joyce", "payed", "brays", "malts", "dregs", "leaky", "runny", "fecal", "woven", "hurls", "jorge", "henna", "dolby", "booty", "brett", "dykes", "rural", "fight", "feels", "flogs", "brunt", "preen", "elvis", "dopey", "gripe", "garry", "gamma", "fling", "space", "mange", "storm", "arron", "hairs", "rogue", "repel", "elgar", "ruddy", "cross", "medan", "loses", "howdy", "foams", "piker", "halts", "jewel", "avery", "stool", "cruel", "cases", "ruses", "cathy", "harem", "flour", "meted", "faces", "hobos", "charm", "jamar", "cameo", "crape", "hooey", "reefs", "denny", "mitts", "sores", "smoky", "nopes", "sooty", "twirl", "toads", "vader", "julep", "licks", "arias", "wrote", "north", "bunks", "heady", "batch", "snaps", "claws", "fouls", "faded", "beans", "wimps", "idled", "pulse", "goons", "noose", "vowed", "ronda", "rajas", "roast", "allah", "punic", "slows", "hours", "metal", "slier", "meaty", "hanna", "curvy", "mussy", "truth", "troys", "block", "reels", "print", "miffs", "busts", "bytes", "cream", "otter", "grads", "siren", "kilos", "dross", "batty", "debts", "sully", "bares", "baggy", "hippy", "berth", "gorky", "argon", "wacko", "harry", "smoke", "fails", "perms", "score", "steps", "unity", "couch", "kelly", "rumps", "fines", "mouth", "broth", "knows", "becky", "quits", "lauri", "trust", "grows", "logos", "apter", "burrs", "zincs", "buyer", "bayer", "moose", "overt", "croon", "ousts", "lands", "lithe", "poach", "jamel", "waive", "wiser", "surly", "works", "paine", "medal", "glads", "gybed", "paint", "lorre", "meant", "smugs", "bryon", "jinni", "sever", "viols", "flubs", "melts", "heads", "peals", "aiken", "named", "teary", "yalta", "styes", "heist", "bongs", "slops", "pouts", "grape", "belie", "cloak", "rocks", "scone", "lydia", "goofs", "rents", "drive", "crony", "orlon", "narks", "plays", "blips", "pence", "march", "alger", "baste", "acorn", "billy", "croce", "boone", "aaron", "slobs", "idyls", "irwin", "elves", "stoat", "doing", "globe", "verve", "icons", "trial", "olsen", "pecks", "there", "blame", "tilde", "milky", "sells", "tangy", "wrack", "fills", "lofty", "truce", "quark", "delia", "stowe", "marty", "overs", "putty", "coral", "swine", "stats", "swags", "weans", "spout", "bulky", "farsi", "brest", "gleam", "beaks", "coons", "hater", "peony", "huffy", "exert", "clips", "riven", "payer", "doped", "salas", "meyer", "dryad", "thuds", "tilts", "quilt", "jetty", "brood", "gulch", "corps", "tunic", "hubby", "slang", "wreck", "purrs", "punch", "drags", "chide", "sulks", "tints", "huger", "roped", "dopes", "booby", "rosin", "outer", "gusto", "tents", "elude", "brows", "lease", "ceres", "laxer", "worth", "necks", "races", "corey", "trait", "stuns", "soles", "teems", "scrip", "privy", "sight", "minor", "alisa", "stray", "spank", "cress", "nukes", "rises", "gusts", "aurae", "karma", "icing", "prose", "biked", "grand", "grasp", "skein", "shaky", "clump", "rummy", "stock", "twain", "zoned", "offed", "ghats", "mover", "randy", "vault", "craws", "thees", "salem", "downy", "sangs", "chore", "cited", "grave", "spinx", "erica", "raspy", "dying", "skips", "clerk", "paste", "moved", "rooks", "intel", "moses", "avers", "staid", "yawls", "blast", "lyres", "monks", "gaits", "floor", "saner", "waver", "assam", "infer", "wands", "bunch", "dryly", "weedy", "honey", "baths", "leach", "shorn", "shows", "dream", "value", "dooms", "spiro", "raped", "shook", "stead", "moran", "ditto", "loots", "tapir", "looms", "clove", "stops", "pinks", "soppy", "ripen", "wench", "shone", "bauds", "doric", "leans", "nadia", "cries", "camus", "boozy", "maris", "fools", "morns", "bides", "greek", "gauss", "roget", "lamar", "hazes", "beefy", "dupes", "refed", "felts", "larry", "guile", "ables", "wants", "warns", "toils", "bathe", "edger", "paced", "rinks", "shoos", "erich", "whore", "tiger", "jumpy", "lamas", "stack", "among", "punts", "scalp", "alloy", "solon", "quite", "comas", "whole", "parse", "tries", "reeve", "tiled", "deena", "roomy", "rodin", "aster", "twice", "musts", "globs", "parch", "drawn", "filch", "bonds", "tells", "droop", "janis", "holds", "scant", "lopes", "based", "keven", "whiny", "aspic", "gains", "franz", "jerri", "steel", "rowel", "vends", "yelps", "begin", "logic", "tress", "sunni", "going", "barge", "blood", "burns", "basks", "waifs", "bones", "skill", "hewer", "burly", "clime", "eking", "withs", "capek", "berta", "cheap", "films", "scoot", "tweed", "sizer", "wheat", "acton", "flung", "ponds", "tracy", "fiver", "berra", "roger", "mutes", "burke", "miked", "valve", "whisk", "runes", "parry", "toots", "japes", "roars", "rough", "irons", "romeo", "cages", "reeks", "cigar", "saiph", "dully", "hangs", "chops", "rolls", "prick", "acuff", "spent", "sulla", "train", "swell", "frets", "names", "anita", "crazy", "sixth", "blunt", "fewer", "large", "brand", "slick", "spitz", "rears", "ogres", "toffy", "yolks", "flock", "gnawn", "eries", "blink", "skier", "feted", "tones", "snail", "ether", "barbs", "noses", "hears", "upset", "awash", "cloud", "trunk", "degas", "dungs", "rated", "shall", "yeahs", "coven", "sands", "susan", "fable", "gunny", "began", "serfs", "balls", "dinky", "madge", "prong", "spilt", "lilly", "brawn", "comet", "spins", "raids", "dries", "sorts", "makes", "mason", "mayra", "royce", "stout", "mealy", "pagan", "nasal", "folds", "libby", "coups", "photo", "mosey", "amens", "speed", "lords", "board", "fetal", "lagos", "scope", "raked", "bonus", "mutts", "willy", "sport", "bingo", "thant", "araby", "bette", "rebel", "gases", "small", "humus", "grosz", "beset", "slays", "steve", "scrap", "blahs", "south", "pride", "heels", "tubes", "beady", "lacey", "genus", "mauls", "vying", "spice", "sexes", "ester", "drams", "today", "comae", "under", "jests", "direr", "yoked", "tempi", "early", "boats", "jesus", "warts", "guppy", "gilda", "quota", "token", "edwin", "ringo", "gaped", "lemon", "hurst", "manor", "arrow", "mists", "prize", "silas", "blobs", "diets", "ervin", "stony", "buddy", "bates", "rabid", "ducat", "ewing", "jaunt", "beads", "doyen", "blush", "thoth", "tiles", "piper", "short", "peron", "alley", "decks", "shunt", "whirs", "cushy", "roils", "betty", "plugs", "woken", "jibes", "foray", "merak", "ruing", "becks", "whale", "shoot", "dwelt", "spawn", "fairs", "dozed", "celts", "blond", "tikes", "sabin", "feint", "vamps", "cokes", "willa", "slues", "bills", "force", "curst", "yokel", "surer", "miler", "fices", "arced", "douse", "hilly", "lucio", "tongs", "togas", "minty", "sagas", "pates", "welsh", "bruno", "decal", "elate", "linux", "gyros", "pryor", "mousy", "pains", "shake", "spica", "pupal", "probe", "mount", "shirk", "purus", "kilns", "rests", "graze", "hague", "spuds", "sweep", "momma", "burch", "maces", "samar", "brace", "riser", "booms", "build", "camel", "flyer", "synge", "sauna", "tonga", "tings", "promo", "hides", "clair", "elisa", "bower", "reins", "diann", "lubed", "nulls", "picks", "laban", "milch", "buber", "stomp", "bosom", "lying", "haled", "avert", "wries", "macon", "skids", "fumed", "ogles", "clods", "antic", "nosey", "crimp", "purge", "mommy", "cased", "taxes", "covet", "clack", "butch", "panty", "lents", "machs", "exude", "tooth", "adore", "shuck", "asses", "after", "terra", "dices", "aryan", "regor", "romes", "stile", "cairo", "maura", "flail", "eaves", "estes", "sousa", "visas", "baron", "civet", "kitty", "freed", "ralph", "tango", "gawks", "cheat", "study", "fancy", "fiber", "musks", "souse", "brims", "claim", "bikes", "venue", "sired", "thymi", "rivas", "skimp", "pleas", "woman", "gimpy", "cawed", "minos", "pints", "knock", "poked", "bowen", "risky", "towel", "oinks", "linus", "heals", "pears", "codas", "inner", "pitch", "harpy", "niger", "madly", "bumpy", "stair", "files", "nobel", "celli", "spars", "jades", "balmy", "kooky", "plums", "trues", "gloss", "trims", "daunt", "tubby", "dared", "wadis", "smell", "darby", "stink", "drill", "dover", "ruler", "laden", "dikes", "layla", "fells", "maker", "joked", "horns", "these", "baize", "spahn", "whens", "edged", "mushy", "plume", "tucks", "spurs", "husky", "dried", "bigot", "pupas", "drily", "aware", "hagar", "newly", "knots", "pratt", "feces", "sabik", "watts", "cooke", "riles", "seamy", "fleas", "dusts", "barfs", "roans", "pawns", "vivid", "kirks", "tania", "feral", "tubae", "horne", "aries", "brits", "combs", "chunk", "stork", "waned", "texan", "elide", "glens", "emery", "autos", "trams", "dosed", "cheri", "baits", "jacks", "whose", "fazed", "matte", "swans", "maxed", "write", "spays", "orion", "traci", "horse", "stars", "strut", "goods", "verge", "scuff", "award", "dives", "wires", "burnt", "dimly", "sleds", "mayan", "biped", "quirk", "sofia", "slabs", "waste", "robby", "mayor", "fatty", "items", "bowel", "mires", "swarm", "route", "swash", "sooth", "paved", "steak", "upend", "sough", "throw", "perts", "stave", "carry", "burgs", "hilts", "plane", "toady", "nadir", "stick", "foist", "gnarl", "spain", "enter", "sises", "story", "scarf", "ryder", "glums", "nappy", "sixes", "honed", "marcy", "offer", "kneel", "leeds", "lites", "voter", "vince", "bursa", "heave", "roses", "trees", "argos", "leann", "grimy", "zelma", "crick", "tract", "flips", "folks", "smote", "brier", "moore", "goose", "baden", "riled", "looks", "sober", "tusks", "house", "acmes", "lubes", "chows", "neath", "vivas", "defer", "allay", "casey", "kmart", "pests", "proms", "eying", "cider", "leave", "shush", "shots", "karla", "scorn", "gifts", "sneer", "mercy", "copes", "faxed", "spurt", "monet", "awoke", "rocky", "share", "gores", "drawl", "tears", "mooed", "nones", "wined", "wrens", "modem", "beria", "hovel", "retch", "mates", "hands", "stymy", "peace", "carat", "coots", "hotel", "karen", "hayed", "mamet", "cuing", "paper", "rages", "suave", "reuse", "auden", "costs", "loner", "rapes", "hazel", "rites", "brent", "pumps", "dutch", "puffs", "noons", "grams", "teats", "cease", "honda", "pricy", "forgo", "fleck", "hired", "silos", "merge", "rafts", "halon", "larks", "deere", "jello", "cunts", "sifts", "boner", "morin", "mimes", "bungs", "marie", "harts", "snobs", "sonic", "hippo", "comes", "crops", "mango", "wrung", "garbs", "natty", "cents", "fitch", "moldy", "adams", "sorta", "coeds", "gilds", "kiddy", "nervy", "slurp", "ramon", "fuzed", "hiker", "winks", "vanes", "goody", "hawks", "crowd", "bract", "marla", "limbs", "solve", "gloom", "sloop", "eaton", "memos", "tames", "heirs", "berms", "wanes", "faint", "numbs", "holes", "grubs", "rakes", "waist", "miser", "stays", "antis", "marsh", "skyed", "payne", "champ", "jimmy", "clues", "fatal", "shoed", "freon", "lopez", "snowy", "loins", "stale", "thank", "reads", "isles", "grill", "align", "saxes", "rubin", "rigel", "walls", "beers", "wispy", "topic", "alden", "anton", "ducts", "david", "duets", "fries", "oiled", "waken", "allot", "swats", "woozy", "tuxes", "inter", "dunne", "known", "axles", "graph", "bumps", "jerry", "hitch", "crews", "lucia", "banal", "grope", "valid", "meres", "thick", "lofts", "chaff", "taker", "glues", "snubs", "trawl", "keels", "liker", "stand", "harps", "casks", "nelly", "debby", "panes", "dumps", "norma", "racks", "scams", "forte", "dwell", "dudes", "hypos", "sissy", "swamp", "faust", "slake", "maven", "lowed", "lilts", "bobby", "gorey", "swear", "nests", "marci", "palsy", "siege", "oozes", "rates", "stunt", "herod", "wilma", "other", "girts", "conic", "goner", "peppy", "class", "sized", "games", "snell", "newsy", "amend", "solis", "duane", "troop", "linda", "tails", "woofs", "scuds", "shies", "patti", "stunk", "acres", "tevet", "allen", "carpi", "meets", "trend", "salty", "galls", "crept", "toner", "panda", "cohen", "chase", "james", "bravo", "styed", "coals", "oates", "swami", "staph", "frisk", "cares", "cords", "stems", "razed", "since", "mopes", "rices", "junes", "raged", "liter", "manes", "rearm", "naive", "tyree", "medic", "laded", "pearl", "inset", "graft", "chair", "votes", "saver", "cains", "knobs", "gamay", "hunch", "crags", "olson", "teams", "surge", "wests", "boney", "limos", "ploys", "algae", "gaols", "caked", "molts", "glops", "tarot", "wheal", "cysts", "husks", "vaunt", "beaus", "fauns", "jeers", "mitty", "stuff", "shape", "sears", "buffy", "maced", "fazes", "vegas", "stamp", "borer", "gaged", "shade", "finds", "frock", "plods", "skied", "stump", "ripes", "chick", "cones", "fixed", "coled", "rodeo", "basil", "dazes", "sting", "surfs", "mindy", "creak", "swung", "cadge", "franc", "seven", "sices", "weest", "unite", "codex", "trick", "fusty", "plaid", "hills", "truck", "spiel", "sleek", "anons", "pupae", "chiba", "hoops", "trash", "noted", "boris", "dough", "shirt", "cowls", "seine", "spool", "miens", "yummy", "grade", "proxy", "hopes", "girth", "deter", "dowry", "aorta", "paean", "corms", "giant", "shank", "where", "means", "years", "vegan", "derek", "tales" };
	unordered_set<string> dict(theSet, theSet + sizeof(theSet)/sizeof(string));
	/*vector<string> res = solution.wordBreak_2(s, dict);
	for (int i = 0; i < res.size(); i++) {
		for (int j = 0; j < res[i].size(); j++)
			cout << res[i] << '|';
		cout << endl;
	}*/
	cout << solution.ladderLength(start, end, dict) << endl;
	

	cin>>stop;
	return 0;
}