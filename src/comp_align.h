#include<bits/stdc++.h>
#include<unistd.h>
#include"tools.h"

void read_align(int k);
void* cal_pre(void* arg);
int dsu(int id1,int id2);
void bas_align(int k);
void* cal_nxt(void* arg);
//read_align and absket align
void align()
{
    //read_align
    pre.resize(rcnt0);
    basket.resize(rcnt0);
    basket_rev.resize(rcnt0);
    //first_rev.resize(rcnt0);
    vecr.clear();
    for(int i=0;i<rcnt0;i++){
        if(isrepeat[i]==0&&hasn[i]==0){
            pre[i]=i,basket[i].pb(i),basket_rev[i].pb(0);
            vecr.pb(i);
        }
        else pre[i]=-1e9;
    }
    
    assert(vecr.size()==rcnt&&rcnt==rcnt0-repeatcnt-ncnt);
    #if testflag
    read_align(4);
    #else
    //read_align(29);
    read_align(27);
    #endif
    cout<<"read_align finish"<<'\n';
    print_time();

    //basket align
    vecr.clear();
    bas_edge.resize(rcnt0);
    first_rev.assign(rcnt0,0);
    for(int i=0;i<rcnt0;i++){//the first edge of every read is outedge(bas_edge[i][0]),others are inedge
        if(isrepeat[i]==0&&hasn[i]==0&&find(i)==i){
            vecr.pb(i);
            bas_edge[i].resize(1);
            bas_edge[i][0]=Edge(i);
        }
    }
    #if testflag
    bas_align(4);
    #else
    //bas_align(29);
    bas_align(27);
    #endif

    //add reverse edge
    //map<int,int> mtmp;
    for(int i=0;i<rcnt0;i++){
        if(isrepeat[i]==1||hasn[i]==1||find(i)!=i) continue;
        if(bas_edge[i][0].bias<=max_bias){
            Edge& tmp=bas_edge[i][0];
            bas_edge[tmp.to].pb(Edge(i,-tmp.bias,tmp.rev));
            first_rev[tmp.to]=first_rev[i]^tmp.rev;

            //mtmp[tmp.bias]++;
        }
    }
    // cout<<"mtmp:\n";
    // for(auto it:mtmp){
    //     cout<<it.first<<' '<<it.second<<'\n';
    // }
    // cout<<'\n';
    cout<<"align finish\n";
    print_time();
}


void read_align(int k)
{
    cout<<"read_align start\n";
    cal_k_minimizer(k);

    sort(vecr.begin(),vecr.end(),[](int id1,int id2){
        if(val[id1]!=val[id2]) return val[id1]<val[id2];
        else return cal_len(id1)>cal_len(id2);
    });
    int tmp;
    mutex_init();
    for(int i=0;i<thread_num;i++){
        if((tmp=pthread_create(&thread_id[i],NULL,cal_pre,&thread_tmp[i]))!=0){
            cout<<"pthread_create ERROR\n";
            exit(0);
        }
    }
    k_minimizer(&thread_tmp[thread_num]);
    for(int i=0;i<thread_num;i++){
        //cout<<i<<'\n';
        void* tmp;
        pthread_join(thread_id[i],&tmp);
    }
}

void* cal_pre(void* arg)
{
    int cnt=0;
    int thread_id=*(int*)arg;
    while(1){
        int curid;
        pthread_mutex_lock(&my_mutex);
        curid=my_id++;
        pthread_mutex_unlock(&my_mutex);
        cnt++;
        int pos=vecr[curid];
        //cout<<"cal_pre: "<<curid<<' '<<pos<<'\n';
        if(curid>=vecr.size()) break;
        for(int i=curid-1;i>=0&&i>=curid-300;i--){
            int ppos=vecr[i];
            if(val[pos]!=val[ppos]) break;
            if(cal_len(pos)!=cal_len(ppos)) break;
            //printf("+++ %d %d %d\n",pos,ppos,check(pos,ppos));
            if(check(pos,ppos,threshold1)){
                //merge two baskets,need lock, the total complexity is O(nlog^2n)
                pthread_mutex_lock(&array_mutex);
                int id1=find(pos),id2=find(ppos);
                pre[id1]=pre[id2]=dsu(id1,id2);
                //first_rev[pos]=iskmersymm[pos];//if multi,then the first_rev is cur rev state
                pthread_mutex_unlock(&array_mutex);
                
                break;
            }
        }
    }
    //printf("thread_id=%d,cnt=%d\n",thread_id,cnt);
    return (void*)0;
}

int dsu(int id1,int id2)
{
    assert(pre[id1]==id1);
    assert(pre[id2]==id2);
    if(id1==id2) return id1;
    if(basket[id1].size()<basket[id2].size()) swap(id1,id2);
    pre[id2]=id1;
    int tmprev=iskmersymm[id1]^iskmersymm[id2];
    for(int i=0;i<basket[id2].size();i++){
        basket[id1].pb(basket[id2][i]);
        basket_rev[id1].pb(basket_rev[id2][i]^tmprev);
    }
    basket[id2].clear();
    return id1;
}


void bas_align(int k)
{
    cout<<"bas_align start,k="<<k<<'\n';
    cal_k_minimizer(k);


    sort(vecr.begin(),vecr.end(),[](int id1,int id2){
        if(val[id1]!=val[id2]) return val[id1]<val[id2];
        else return cal_len(id1)>cal_len(id2);
    });
    
    nxt_bas.assign(vecr.size(),1e9);
    int tpre=1e9;
    for(int i=vecr.size()-2;i>=0;i--){
        if(val[vecr[i]]!=val[vecr[i+1]]) tpre=i+1;
        else if(cal_len(vecr[i])!=cal_len(vecr[i+1])) tpre=i+1;
        nxt_bas[i]=tpre;
    }

    #if testflag
    cout<<"vecr:\n";
    for(int i=0;i<vecr.size();i++) cout<<vecr[i]<<' ';
    cout<<'\n';
    #endif
    int tmp;
    
    mutex_init();

    for(int i=0;i<thread_num;i++){
        if((tmp=pthread_create(&thread_id[i],NULL,cal_nxt,&thread_tmp[i]))!=0){
            cout<<"pthread_create ERROR\n";
            exit(0);
        }
    }
    cal_nxt(&thread_tmp[thread_num]);
    for(int i=0;i<thread_num;i++){
        void* tmp;
        pthread_join(thread_id[i],&tmp);
    }
    
    cout<<"bas_align finish\n";
}

void* cal_nxt(void* arg)
{
    int id=*(int*) arg;
    while(1){
        int curid;
        pthread_mutex_lock(&my_mutex);
        curid=my_id++;
        pthread_mutex_unlock(&my_mutex);
        //printf("pthread_id=%d,curid=%d\n",id,curid);
        if(curid>=vecr.size()) break;
        int pos=vecr[curid];
        //printf("cal_nxt: %d %d %d\n",curid,pos,nxt_bas[curid]);
        for(int i=nxt_bas[curid];i<vecr.size()&&i<=nxt_bas[curid]+300;i++){
            int npos=vecr[i];
            if(val[pos]!=val[npos]) break;
            if(cal_len(pos)-cal_len(npos)>max_bias) break;
            assert(cal_len(pos)>cal_len(npos));
            if(check(pos,npos,threshold2)){
                //merge two baskets,need lock, the total complexity is O(nlog^2n)
                pthread_mutex_lock(&array_mutex);
                Edge tmpe=Edge(npos,cal_len(pos)-cal_len(npos),iskmersymm[pos]^iskmersymm[npos]);
                int tmp=bas_edge[pos][0].update(tmpe);
                if(tmp) first_rev[pos]=iskmersymm[pos];
                //else ,bias too large ,discard
                pthread_mutex_unlock(&array_mutex);//can't merge,but this can occur only few cases
                break;
                
            }
        }
    }
    return (void*)0;
}