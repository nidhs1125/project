#include<bits/stdc++.h>
#include"tools.h"

void cal_repeat();
void cal_n();

void prework(ifstream& fin)
{
    hash_base1.assign(max_str_length,1);
    hash_base2.assign(max_str_length,1);
    invhash_base1.assign(max_str_length,1);
    invhash_base2.assign(max_str_length,1);
    for(int i=1;i<max_str_length;i++){
        hash_base1[i]=hash_base1[i-1]*step1%mod1;
        invhash_base1[i]=inv(hash_base1[i],mod1);
        hash_base2[i]=hash_base2[i-1]*step2%mod2;
        invhash_base2[i]=inv(hash_base2[i],mod2);
    }

    string tmp;
    int cnt=0;
    vecr.clear();
    while(getline(fin,tmp)){
        if(tmp[0]=='$') break;
        cnt++;
        if(cnt%4==2){
            vread.pb(tmp);
            vecr.pb(rcnt++);
        }
    }
    cout<<rcnt<<"reads in total\n";
    rcnt0=rcnt;
    read_len=vread[0].length();
    maxdiscnt=1.0*read_len/threshold_ratio_re;
    for(int i=1;i<rcnt0;i++) assert(vread[i].length()==read_len);
    
    if(read_len>max_str_length){
        printf("read_len >max_str_length, terminated\n");
        exit(0);
    }
    print_time();
    //exclude the same reads
    isrepeat.assign(rcnt0,0);
    issymmrepeat.assign(rcnt0,0);
    repeatid.assign(rcnt0,0);
    cal_repeat();//find repeat read
    cout<<"cal_repeat finish\n";
    hasn.assign(rcnt0,0);
    cal_n();
    cout<<"prework finish\n";
}

void cal_repeat()
{
    cout<<"start cal_repeat\n";
    vector<int> vec;
    hash_val.resize(rcnt0);
    for(int i=0;i<rcnt0;i++){
        hash_val[i]=get_hash_val(vread[i]);
        vec.pb(i);
    }
    sort(vec.begin(),vec.end(),[](int& s1,int& s2){
        if(hash_val[s1].first!=hash_val[s2].first) return hash_val[s1].first<hash_val[s2].first;
        else return hash_val[s1].second<hash_val[s2].second;
    });
    for(int i=0;i<rcnt0;i++){
        if(isrepeat[vec[i]]==1) continue;
        //try to find repeat read
        int j=i+1;
        while(j<rcnt0&&hash_val[vec[i]]==hash_val[vec[j]]){
            if(isrepeat[vec[j]]!=0) {j++;continue;}
            repeatcnt++;
            isrepeat[vec[j]]=1;
            issymmrepeat[vec[j]]=0;
            repeatid[vec[j]]=vec[i];
            j++;
        }
        //try to find symmetric repeat read
        string tmps=cal_symm(vread[vec[i]]);
        puu tmpp=get_hash_val(tmps);
        //binary search
        int l=0,r=rcnt0-1,d=0;
        while(l<r){
            d=l+r>>1;
            if(tmpp.first==hash_val[vec[d]].first){
                if(tmpp.second<=hash_val[vec[d]].second) r=d;
                else l=d+1;
            }
            else if(tmpp.first<hash_val[vec[d]].first) r=d;
            else if(tmpp.first>hash_val[vec[d]].first) l=d+1;
        }
        j=l;
        while(j<rcnt0&&tmpp==hash_val[vec[j]]){
            if(isrepeat[vec[j]]!=0) {j++;continue;}
            assert(isrepeat[vec[j]]==0);
            repeatcnt++;
            isrepeat[vec[j]]=1;
            issymmrepeat[vec[j]]=1;
            repeatid[vec[j]]=vec[i];
            j++;
        }
        if(isrepeat[vec[i]]==1) isrepeat[vec[i]]=0,repeatcnt--;//reverse symm with itself
    }
    vecr.clear();
    for(int i=0;i<rcnt0;i++) if(isrepeat[i]==0) vecr.pb(i);
    //sort(vecr.begin(),vecr.end(),cmp3);
    rcnt-=repeatcnt;
    assert(vecr.size()==rcnt);
}

//find all read(not repeat) with n.
void cal_n()
{
    cout<<"start cal_n\n";
    sort(vecr.begin(),vecr.begin()+rcnt,[](int& r1,int& r2){
        return vread[r1].find('N')<vread[r2].find('N');//if don't has n,result=-1,in first
    });
    for(int i=rcnt-1;i>=0;i--){
        if(vread[vecr[i]].find('N')!=string::npos){
            rcnt--,ncnt++;
            hasn[vecr[i]]=1;
        }
    }
    vecr.resize(rcnt);

}