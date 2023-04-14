#include<bits/stdc++.h>
#include<unistd.h>
#include <sys/time.h>
using namespace std;
#define pb push_back
vector<string> vec1,vec2;
char symm(char c,bool is)
{
    char ret;
    if(is==0){
        ret=c;
    }
    else{
        switch(c){
            case 'A':
                ret='T';
                break;
            case 'T':
                ret='A';
                break;
            case 'C':
                ret='G';
                break;
            case 'G':
                ret='C';
                break;
            default:
                ret='N';
        }
    }
    
    return ret;
}
string cal_symm(string& s)
{
    string tmp="";
    for(int p=s.length()-1;p>=0;p--){
        tmp+=symm(s[p],1);
    }
    return tmp;
}
signed main()
{
    while(1){
        string s;
        cin>>s;
        cout<<cal_symm(s);
    }
    
    return 0;
}