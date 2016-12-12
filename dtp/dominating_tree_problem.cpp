#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<iostream>
#include<algorithm>
using namespace std;
#define MAX 400
#define MAXCOST 0x7fffffff
double graph[MAX][MAX];
double t[MAX][MAX];//备份于操作
double tree[MAX][MAX];//树形结构
double tpp[MAX][MAX];
double Wtmp[MAX][MAX];
double tmp[MAX][MAX];
double treetmp[MAX][MAX]; 
int v[MAX];//点集
int d[MAX];//支配集序列
int a[MAX];//real premution
int b[MAX];//存储的一系列顶点的排序b[]
int r[MAX];
int x[MAX];
int xsharp[MAX];
int xsharp2[MAX];
int xpp[MAX];
int xstar[MAX];
int len_maxvnum;
int lenrm_maxvnum;
int len_jl;
int vnum;
int len;
int len_treetmp;//len(treetmp)
int lenprim;//a的长度
int lenprim2;//a的长度
int lenshak;//剪枝后的x的长度
int lenshak2;
int leni;//b的长度
double f,fpp,fstar;
double p=0.3;
int kk,kmax=20,kmin=2;
int m,n;
int len_x;
///v_len在变，应该用一个const变量
///这里的lenth都有问题，传的应该是最大顶点数，而不是一共的顶点的个数
void MakePermutation()
{
    int flag[m];
    memset(flag,0,sizeof(flag));
    for(int i=0;i<m;i++)
        for(int j=0;j<m;j++)
        {
            if(graph[i][j]!=0)
            {
                flag[i]=flag[j]=1;
            }
        }
    int count=0;
    for (int i = 0; i < m; ++i)
    {
        
        if(flag[i]==1)
            x[count++]=i;
    }
    len=count;//初始化的x的长度
}
//是不是应该是最大点而不是一共多少个点,posit 的位置！！！
void shaking(int px[],int k)//k=8的时候不行啊
{
    
    cout<<"------shaking-----"<<endl;
    cout<<"------len_x-------"<<endl;
    cout<<"len_x = "<<len_x<<endl;
    cout<<"最大顶点是lenrm_maxvnum = "<<lenrm_maxvnum<<endl;
    cout<<"------px[i]-------"<<endl;
    for (int i = 0; i < m; ++i)
    {   cout<<px[i]<<" ";
        xpp[i]=px[i];
    }
    cout<<endl;
    int tmp[MAX];
    int posit[MAX];
    memset(posit,0,sizeof(posit));
    memset(tmp,-1,sizeof(tmp));
    int i=0;
    ///传进来的是选k个顶点，比如k=3，选择点的位置1，3，6，
    ///然后就把1位置的顶点随机放在3，6位置，
    ///3位置的顶点随机放在1，6位置
    //tmp[0]=rand()%len_x

    while(k)
    {
        //cout<<"k = "<<k<<endl;
        //随机在x里选取k个点，记录
        int flgg=0;
        int t = rand()%m;
        int tt=xpp[t];
        cout<<"i = "<< i <<endl;
        for (int j = 0; j < i; ++j)
        {
            if(tmp[j]==tt)
            {
                flgg=1;
                break;
            }
        }
        if(flgg==0)
        {
        tmp[i++]=xpp[t];//x里面存的是点,
        posit[tt]=1;
        k--;
        }
    }
    cout<<"随机选出的数 "<<endl;
    for (int j = 0; j < i ; ++j)
    {
       cout<<tmp[j]<<" ";
    }
    cout<<endl;
    ///传进来的是选k个顶点，比如k=3，选择点的位置1，3，6，
    ///然后就把1位置的顶点随机放在3，6位置，
    ///3位置的顶点随机放在1，6位置
    
    int tswap;
    int tf;
    int tn;
    int c;//不断更新tmpp
    int  xppt[MAX];
    int tmpp[MAX];
    memset(tmpp,-1,sizeof(tmpp));
    memset(xppt,-1,sizeof(xppt));
    for (int j = 0; j < i; ++j)
    {
       tmpp[j]=tmp[j];
    }
    c=i;
    //cout<<"c = "<<rand()%c<<endl;
    int zz;

   for (int j = 0; j < i ; ++j)
   {
       //cout<<"j = "<<j<<endl;
       if(c==0)break;
       if(c==1)
            {
                for (int z = 0; z <i; ++z)
                {
                  if(tmpp[z]!=-1)
                    {
                        tn=tmpp[z];//得到选出的 位置
                        zz=z;
                        c--;//这里需不需要一个break；
                    }
                    
                }
                  
            }
       while(c>1)
       {
            int count =-1;
            //cout<<"c = "<<c<<endl;  
            tf=rand()%c;//c个点中随机选取一个点
        for (int z = 0; z <i; ++z)
        {
           if(tmpp[z]!=-1)//如果某个位置存在
            count++;//操作加一
            if(count==tf)//如果位置＝到了tf就说明找到了随机的位置
                {
                    tn=tmpp[z];//得到选出的 位置
                    zz=z;
                    //位置个数减
                    break;
                }
        }
        
        //tn=tmpp[tf];
        //cout<<" tn = "<<tn;
        if(tmp[j]!=tn&&posit[tn]!=0)//如果要摆放的点和选出的位置不同并且tn的位置存在
            {
                tmpp[zz]=-1;//减去这个位置
                posit[tn]=0;
                    c--;
                break;//选出位置
            }    
       }

       for (int z = 0; z < m; ++z)//lenrm_maxvnum
       {
          if(xpp[z]==tn)//找到位置
          {
            xppt[z]=tmp[j];//赋值
            posit[tn]=0;
          }
       }
   }
   
    for (int j = 0; j < m; ++j)
    {
        if(xppt[j]!=-1)
            xpp[j]=xppt[j];
    }
    cout<<"shaking后的数组 "<<endl;
    for (int j = 0; j < m; ++j)
    {
       cout<<xpp[j]<<"_";
    }
           cout<<endl;
    //这里用不用做个x的替身
           cout<<"------shaking_end-----"<<endl;
}
double wsum(int lenth,double g[MAX][MAX])//lenth传入最大顶点＋1
{
    double sum=0;
    for (int i = 0; i < lenth; ++i)
    {
        for (int j = 0; j < lenth; ++j)
        {
            if(g[i][j]!=0)
                sum+=g[i][j];
        }
    }
    return sum/2;

}
//传入的是最大顶点
double removelef(int lenth ,double g[MAX][MAX])//剪枝后最小生成树顺序存在x数组里，tree存的是生成树，都需要备份xsharp
{
    f=0;//有个数组需要初始化
    cout<<"----------rl---------"<<endl;
    memset(r,0,sizeof(r));
    for (int i = 0; i < lenth; ++i)
    {
        for (int j = i+1; j < lenth; ++j)
        {
            if(d[i]!=0&&d[j]!=0)
            {
                int t1=i;
                int t2=j;//a[j]不是顺序的？？？这要改
                //cout<<"t1="<<t1<<" t2="<<t2<<endl;
                if(g[t1][t2]!=0)
                {
            //cout<<"t1="<<t1<<" t2="<<t2<<endl;
            //cout<<"g="<<g[t1][t2]<<endl;
                    r[t1]++;
                    r[t2]++;
                }
            }
            
        }
    }
    lenshak=0;
    
    for (int i = 0; i < lenth; ++i)
    {
        if(r[i]<=1)
        {

            for (int j = 0; j < lenth; ++j)//这个界是lenprim？
            {
              g[i][j]=g[j][i]=0; 
            }
            r[i]=0;
        }
        
       if(r[i]>1)
        {
            xsharp[lenshak++]=i;
            //cout<<"d"<<i<<endl;
        }

    }
    lenrm_maxvnum=xsharp[lenshak-1];
    //cout<<"lenrm_maxvnum = "<<len_maxvnum<<endl;
    for (int i = 0; i < lenth; ++i)
    {
        for (int j = 0; j < lenth; ++j)
        {
           //if(g[i][j]!=0)
           //{
            //cout<<" i= "<<i<<" j= "<<j;
           //}
            f+=g[i][j];
        }
        //cout<<endl;
    }
    cout<<" remove leaf f = "<<f<<endl;
    return f;
}
void prim1(int b[],int lenth,double g[MAX][MAX])//这里的lenth估计有问题
{
    cout<<"----------prim1---------"<<endl;
    int mst[MAX];
    double lowcost[MAX];//equal 0?join:not
    double sum =0;
    for (int i = 0; i < MAX; ++i)
    {
        lowcost[i]=MAXCOST;
    }
    //init
    memset(d,0,sizeof(d));
    memset(a,0,sizeof(a));
    memset(mst,-1,sizeof(mst));
    for (int i = 0; i < MAX; ++i)
    {
       for (int j = 0; j < MAX; ++j)
       {
           t[i][j]=MAXCOST;
           tree[i][j]=0;
       }
    }

    for (int i = 0; i < lenth-1; ++i)//赋值
    {
        for (int j = i+1; j < lenth; ++j)
        {
            int t1=b[i];
            int t2=b[j];
            if(g[t1][t2]!=0)
            {
                t[t1][t2]=g[t1][t2];
                t[t2][t1]=g[t1][t2];
                //cout<<t1<<" "<<t2<<" "<<t[t1][t2]<<endl;
            }
        }
        
    } 
    int fp=b[0];//fitst point
    for (int i = 1; i<lenth; ++i)
    {
        int dp=b[i];
        lowcost[dp]=t[fp][dp];//以dp为终点
        mst[dp]=b[0];//lowcost的起点
    }  
    mst[fp]=0;//把d0加入mst
    for (int i = 1; i <lenth; ++i)
    {
        double min = MAXCOST;
        int minid=-1;
        for (int j = 1; j<lenth; ++j)
        {
            int dp=b[j];
            if(lowcost[dp]<min && lowcost[dp]!=0)
            {
                min =lowcost[dp];
                minid=dp;
            }
        }
        int ttt=mst[minid];
        tree[ttt][minid]=min;
        tree[minid][ttt]=min;
        //构造新的序列？？？这里有问题
        d[ttt]=1;
        d[minid]=1;
        //cout << "V" << mst[minid] << "-V" << minid << "=" << min << endl; 
        sum+=min;
        lowcost[minid]=0;
        for (int j = 1; j<lenth ; ++j)
        {
            int dp=b[j];

            if(t[minid][dp]<lowcost[dp])
            {
                lowcost[dp]=t[minid][dp];
                mst[dp]=minid;

            }
        }
    }
    int count=0;
    for (int i = 0; i < MAX; ++i)
    {
        if(d[i]!=0)
                a[count++]=i;

    }
    lenprim=count;
    len_maxvnum=a[count-1];
    cout<<"len_maxvnum="<<len_maxvnum<<endl;
    
}
int judgeds(int px[],int ppx[],int m,int lenth,double g[MAX][MAX])
{
    //judgeds(b,x,i,len,g)==1
    int fg = 1;
    int flag[lenth];
    for (int i = 0; i < lenth; ++i)
    {
        flag[i]=-1;
    }//m 是b里的顶点个数
    int ttt[MAX];
    memset(ttt,0,sizeof(ttt));
    for (int i = 0; i < m; ++i)
    {
        if(px[i])
        {
            int t=px[i];
            ttt[t]=1;
            flag[t]=1;
        }
            
    }
    for(int i=0;i<lenth;i++)//应该是最大顶点＋1
    {
        

        for (int j = 0; j<m; ++j)//b里面的边
        {
            
            int t1=px[j];//b
            int t =ppx[i];//x
             if(g[t][t1]!=0&&ttt[t]==0)
             {
                flag[i]=1;//外边的顶点被访问到了
             }
        }
         if(flag[i]!=1)
            {
                fg=0;
                break;
            }
    }
    return fg;
}
//连通性
bool Warshall(int b_tmp[],int b[],int mm,int n,double g[MAX][MAX])//vertex num
{
    //memset(Wtmp, MAXCOST, sizeof(Wtmp)); // important

    for (int i = 0; i < MAX; ++i)
    {
        for (int j = 0; j < MAX; ++j)
        {
            Wtmp[i][j]=MAXCOST;
        }
    }
    //输入a数组
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <n; ++j)
        {
            if(b_tmp[i]!=0&&b_tmp[j]!=0)
           {
            if(g[i][j]!=0)
            Wtmp[i][j]=Wtmp[j][i]=g[i][j];
            }
        }
    }
    //初始化tmp数组
    for (int i = 0; i < MAX; ++i)
    {
        for (int j = 0; j < MAX; ++j)
        {
            tmp[i][j]=0;
        }
    }
    //memset(tmp, 0, sizeof(tmp));
    for(int i=0;i<n;i++)//这里的点都穿错了
    {
        for(int j=0;j<n;j++)
        {
            if(Wtmp[i][j] < MAXCOST) // 代表不通
                tmp[i][j] = tmp[j][i] = 1;
            else
                tmp[i][j] = tmp[j][i] = 0;
        }
        tmp[i][i] = 1;
    }
    
 for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(tmp[i][j])
            {
                for(int k=0;k<n;k++)
                {
                    if(tmp[k][i])
                        tmp[k][j] = tmp[j][k] = 1;
                }
            }
        }
    }
    for(int i=0;i<mm;i++)
        for(int j=0;j<mm;j++)
        {
            int t1=b[i];
            int t2=b[j];
            if(!tmp[t1][t2])
                return false;
        }
            
    return true;    
}

//求支配树，这里的问题是传参数。。。能不能用全局变量
void dstree(int le,int px[],double g[MAX][MAX])//传入长度和数组
{
    //dstree(v_len,xpp,graph);
    cout<<"----------dstree---------"<<endl;
    cout<<"-------输出传入的参数px---------"<<endl;
    for (int i = 0; i < le; ++i)
    {
        cout<<px[i]<<" ";
    }
    cout<<endl;
    memset(b,0,sizeof(b));
    int i=0;
    while(1)
    {
        b[i]=px[i];//可以不用全局？
        i++;
        if(judgeds(b,x,i,len,g)==1)
        {
            break;
        }
    }
    cout<<"-------输出judgede b---------"<<endl;
    for (int j = 0; j < i; ++j)
    {
        cout<<b[j]<<" ";
    }
    cout<<endl;
    int b_tmp[MAX];
    memset(b_tmp,0,sizeof(b_tmp));
    while(1)
    {
        int b_count=0;
        for (int j = 0; j < i; ++j)
        {
                int t1=b[j];
                b_tmp[t1]=1;
                if(b[j]>b_count)
                    b_count=b[j];
        }
        //cout<<"b_count = "<<b_count<<endl;
        if(Warshall(b_tmp,b,i,b_count+1,g)==0&&i<le)
        {
            b[i]=px[i];
            i++;
            
        }
        else
        {
            break;
        }
    }
    leni=i; 
    cout<<"-------输出b---------"<<endl;
    for (int j = 0; j < i; ++j)
    {
        cout<<b[j]<<" ";
    }
    cout<<"i = "<<leni<<endl;  
    
}

void localsearch()//int xpp[])
{
    cout<<"＊＊＊＊局部搜索一开始的dstree"<<endl;
    cout<<"xpppppppppp"<<endl;
    for (int i = 0; i < m; ++i)
    {
        cout<<xpp[i]<<" ";
    }
    cout<<endl;
    int  v_len=m;//xpp's lenth
    dstree(v_len,xpp,graph);//t''
    prim1(b,leni,graph);//得到a＝t''和treetmp
    cout<<"leni = "<<leni<<endl;
    cout<<"lenprim = "<<lenprim<<endl;
    fpp=wsum(len_maxvnum+1,tree);//得到xsharp
    int impr =1;
    while(1)
    {
        if(impr==0)
            break;
        int xtmp[MAX];
        memset(xtmp,0,sizeof(xtmp));
        impr=0;
        int flag=0;
        int t_len=lenprim;//tpp's lenth
        for (int i = 0; i < t_len-1; ++i)//lenprim2是当前树的顶点数
        {
            for (int j = v_len-1; j>=t_len; --j)//lenprim时x‘’的顶点数
            {
                for(int z=0;z < v_len;z++)
                {
                    xtmp[z]=xpp[z];//赋值
                }
                xtmp[i]=xpp[j];
                xtmp[j]=xpp[i];
                cout<<"*****在局部搜索里的dstree"<<endl;
                dstree(v_len,xtmp,graph);
                prim1(b,leni,graph);//tree,lenprim都变了
                double ftmp=wsum(len_maxvnum+1,tree);
                cout<<"ftmp = "<<ftmp<<endl;
                cout<<"fpp = "<<fpp<<endl;
                if(ftmp<fpp)
                {
                    cout<<"ftmp<fpp"<<endl;
                    fpp=ftmp;
                    //memset
                    cout<<"xpp[]"<<endl;
                    for(int z=0;z<v_len;z++)
                    {
                        xpp[z]=xtmp[z];
                       // cout<<xpp[z]<<" ";

                    }

                    cout<<endl;
                    for(int i = 0; i < MAX; ++i)
                    {
                     for (int j = 0; j < MAX; ++j)
                        {
                             tpp[i][j]=0;
                         }
                    }
                    len_jl=m;
                    for(int i = 0; i < len_jl+1; ++i)
                    {
                     for (int j = 0; j < len_jl+1; ++j)
                        {
                             tpp[i][j]=tree[i][j];
                         }
                    }
                    cout<<"tpp"<<endl;
                    memset(xstar,-1,sizeof(xstar));
                        for (int i = 0; i < len_jl+1; ++i)
                        {
                        for (int j = 0; j <len_jl+1; ++j)
                        {
                            if(tpp[i][j]!=0)
                                {
                                    xstar[i]++;
                                    xstar[j]++;
                                }
                        }
                        }
                        cout<<" xstar ＝＝＝＝＝";
                        for (int i = 0; i < len_jl+1; ++i)
                        {
                            if(xstar[i]!=-1)
                                cout<< i <<" ";
                        }
                        cout<<endl;
                        int xtest[MAX];
                        int ctest=0;
                        memset(xtest,-1,sizeof(xtest));
                         for (int i = 0; i < len_jl+1; ++i)
                         {
                             if(xstar[i]!=-1)
                                 {
                                        cout<< i <<" ";
                                        xtest[ctest++]=i;
                                    }
                
                            }
        cout<<endl;
        int trues=-1;
        cout<<" true or false"<<endl;
        trues=judgeds(xtest,x,ctest,len,graph);
        cout<<trues<<endl;
                      for (int i = 0; i < len_jl+1; ++i)
                        {
                        for (int j = 0; j <len_jl+1; ++j)
                        {
                            cout<<tpp[i][j]<<" ";

                        }
                        cout<<endl;
                        }
                        

                    t_len=lenshak;
                    impr=1;
                    flag = 1;
                    break;
                }
            }
            if(flag==1)
               {
                cout<<"跳到repeat"<<endl;
                break;

               } 
        }
    }
    cout<<"--------xpp---------"<<endl;
    for (int i = 0; i < v_len; ++i)
    {
       cout<<xpp[i]<<" "<<endl;
    }
    cout<<"local search 得到的解"<<endl;
    cout<<"fpp = "<<fpp<<endl;
}
//这里是猪函数的大框
void mast()
{
    MakePermutation();
    //dstree(len,x,graph);
    cout<<"------初始化－－－－－"<<endl;
    prim1(x,len,graph);//参数问题
    fstar=removelef(lenprim,tree)/2;

    cout<<"fstar="<<fstar<<endl;
    len_x=lenshak;
    
    for (int i = 0; i < len_x; ++i)
    {
       //xstar[i]=xsharp[i];
       xsharp2[i]=xsharp[i];
    }
    int count =len_x;
            for (int i = 0; i < m; ++i)
            {
                 int  flagg =0;
                 for (int j = 0; j < len_x; ++j)
                    {
                    if (x[i]==xsharp[j])
                    {
                        flagg=1;
                        break;
                    }
                    }
        
                    if(flagg ==0)
                    xsharp2[count++]=x[i];
            }
    for (int i = 0; i < m; ++i)
    {
        cout<<xsharp2[i]<<" ";
    }
    cout<<endl;
    //memset(xstar,-1,sizeof(xstar));


        kk=kmin;
        if(kmax>m)
        {
            kmax=m+1;
        }
        while(kk<kmax)
        {
            cout<<"------------k--------"<<endl;
            cout<<"kk="<<kk<<endl;
            shaking(xsharp2,kk);//传入xsharp
            
            localsearch();
            int flag=0;
            if(fpp<fstar)
            {
                
                /*for (int i = 0; i < len_x; ++i)//改
                {
                    xstar[i]=xpp[i];
                }
                memset(xstar,-1,sizeof(xstar));
                        for (int i = 0; i < len_jl; ++i)
                        {
                        for (int j = 0; j <len_jl; ++j)
                        {
                            if(tpp[i][j]!=0)
                                {
                                    xstar[i]++;
                                    xstar[j]++;
                                }
                        }
                        }*/
                cout<<"fstar fuzhi"<<endl;
                /*for (int i = 0; i < len_x; ++i)//改
                {
                    cout<<xpp[i]<<" ";
                }*/
                
                cout<<endl;

                fstar=fpp;
                flag=1;
            }
            else if(fpp==fstar)
            {
                double tmp=rand()%10/(double)10;
                cout<<"随机生成的概率"<<endl;
                cout<<"p_tmp= "<<tmp<<endl;
                if(tmp<p)
                {
                    //memset(xstar,-1,sizeof(xstar));
                    for (int i = 0; i < m ; ++i)
                    {
                    xstar[i]=xpp[i];
                    }
                    /*
                    memset(xstar,-1,sizeof(xstar));
                        for (int i = 0; i < len_jl; ++i)
                        {
                        for (int j = 0; j <len_jl; ++j)
                        {
                            if(tpp[i][j]!=0)
                                {
                                    xstar[i]++;
                                    xstar[j]++;
                                }
                        }
                        }
                        */
                        cout<<"fstar fuzhi2"<<endl;
               
                cout<<endl;
                    //cout<<"fstar fuzhi2"<<endl;
                    fstar=fpp;
                    flag=1;
                }
            }
            if(flag==0)
                kk++;
        }
        cout<<"-----end----"<<endl;
        cout<<"-----end----"<<endl;
        
        cout<<endl;
        cout<<fstar<<endl; 
}
int main()
{
    int i,j,z;
    srand(time(NULL));
    //int m,n;//m个定点n个边
    double weight;
    cin>>m>>n;//init
    for(i=0;i<m;i++)
    {
        for(j=0;j<m;j++)
        {
            graph[i][j]=0;
            t[i][j]=0;
        }
        v[i]=0;
       // b[i]=NULL;
    }//复制
    for (z = 0; z <n; ++z)
    {
        cin>>i>>j>>weight;
        graph[i][j]=weight;
        graph[j][i]=weight;
    }
    mast();
    return 0;
}
