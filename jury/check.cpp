#include<bits/stdc++.h>
#include<unistd.h>
#include <sys/time.h>
using namespace std;
#define pb push_back
vector<string> vec1,vec2;
bool check()
{
    ifstream fin1;
    ifstream fin2;
    fin1.open("test.fastq");
    string tmp;
    int cnt=0;
    while(getline(fin1,tmp)){
        cnt++;
        if(cnt%4==2){
            vec1.pb(tmp);
        }
    }
    fin2.open("decomp.fastq");
    while(getline(fin2,tmp)){
        vec2.pb(tmp);
    }

    if(vec1.size()!=vec2.size()){
        cout<<"No1\n";
        cout<<vec1.size()<<' '<<vec2.size()<<'\n';
        return 0;
    }
    sort(vec1.begin(),vec1.end());
    sort(vec2.begin(),vec2.end());
    for(int i=1;i<vec1.size();i++){
        if(vec1[i]!=vec2[i]){
            cout<<"No2\n";
            return 0;
        }
    }
    return 1;
}
signed main()
{
    int cnt=0;
    while(1){
        //system("./random.out");
        
        //system("./main.out ");
        //system("./main.out -d");
        
        if (!check()){
        	printf("testcase %d : WA\n",cnt);
        	break;
		}
		else{
			printf("testcase %d : AC\n",cnt++);
		}
        break;
    }
    return 0;
}