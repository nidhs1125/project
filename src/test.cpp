
#include<bits/stdc++.h>
#include<unistd.h>
#include<pthread.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
#define mod 998244353ll
ll qow(ll a,ll p) {ll ans=1;for(;p;a=a*a,p>>=1) if(p&1) ans=ans*a;return ans;}
void* get_hash_value(void* arg)
{
    int x=*(int* )arg;
    cout<<"get_hash_value: "<<x<<'\n';
    return ((void*)0);
}
int main(int argc,char** argv)
{
    // string str="AAAAAAAAAAAAAAAAAAAA";
    // string str2="AAAAAAAAAAAAAAAAAAA";
    // cout<<get_hash_value(str)<<' '<<get_hash_value(str2);
    int tmp;
    int x=1;
    pthread_t pid_i;
    cout<<"+++"<<1<<' '<<x<<'\n';
    if((tmp=pthread_create(&pid_i,NULL,get_hash_value,(void*)&x))!=0){
        cout<<"pthread_create ERROR\n";
        exit(0);
    }
    int y=2;
    cout<<"+++"<<2<<' '<<y<<'\n';
    if((tmp=pthread_create(&pid_i,NULL,get_hash_value,(void*)&y))!=0){
        cout<<"pthread_create ERROR\n";
        exit(0);
    }
    sleep(1);
}

/*
#include <iostream>
#include <pthread.h>
#include<unistd.h>
using namespace std;

void* thr_fn(void* arg)
{
        int i = *(int*)arg;
        cout << i << endl;

        return ((void*)0);
}

int main()
{
        pthread_t tid;
        int j = 2;
        pthread_create(&tid, NULL, thr_fn, &j);
        sleep(1);
        return 0;
}
*/