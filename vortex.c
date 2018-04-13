#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define MAX 100		//矩阵大小
#define J 1			//J的值
#define KT 1.5		//KT的值
#define STEP 5			//采样步长
#define TIME 10		//最大的关联步长

struct record{
    int number;
    double avg;	//平均涡旋密度
    double mag;	//平均磁化强度
};

struct exp{
    int k;
    int x;
    int y;
    double p1;
    double p;
    double angle1;
    double angle2;
    double vortex1;
    double vortex2;
};
struct val{
    int x;
    int y;
    double value;
};

void print(double array[MAX+2][MAX+2],double vortex);
void shellsort(struct record s[],int len);
struct exp cal(double array[MAX+2][MAX+2]);
void statisticcount(double array[MAX+2][MAX+2],double *vortex,double *vortexsum,double *mag,double *magsum);
void borders(double array[MAX+2][MAX+2],int x,int y);
double check(double array[MAX+2][MAX+2]);

struct record jilu[5000];


//选择模式//
int main(){
    int type;
    int i,j;
    double array[MAX+2][MAX+2];
    void timetrace(double array[MAX+2][MAX+2]);
    void counttrace(double array[MAX+2][MAX+2]);
    void relation(double array[MAX+2][MAX+2]);
    void point(double array[MAX+2][MAX+2]);
    const double pi=acos(-1.0);
    srand((unsigned)time(NULL));
    for(i=0;i<MAX+2;++i)
        for(j=0;j<MAX+2;++j)
            array[i][j]=0.0;		//有序初态
    //array[i][j]=(rand()/(RAND_MAX+1.0))*pi*2;		//无序初态
    printf("请选择模式：\n（1）随时间演化\n（2）重复多次求总平均\n（3）空间关联函数\n（4）取点计算平均磁化强度\n请选择：");
    scanf("%d",&type);
    if(type==1)
        timetrace(array);
    else if(type==2)
        counttrace(array);
    else if(type==3)
        relation(array);
    else if(type==4)
        point(array);
    return 0;
}

//取点计算//
void point(double array[MAX+2][MAX+2]){
    int n;
    int k=0;
    double vortex,vortexsum=0.0;
    double mag,magsum=0.0;
    printf("输入循环次数(小于%d)：",1000*MAX*MAX*STEP);
    scanf("%d",&n);
    while(k<n){
        cal(array);
        if(k%(MAX*MAX*STEP)==0)
            statisticcount(array,&vortex,&vortexsum,&mag,&magsum);
        k++;
    }
    printf(" 步数               涡旋平均值            磁化平均值            大小顺序\n");
    for(k=0;k<(n/(MAX*MAX*STEP));k++){		//打印随时间演化的值
        printf("   %d          %.16f         %.16f      %3d\n",jilu[k].number,jilu[k].avg,jilu[k].mag,k+1);
    }
}
//随时间的演化//
void timetrace(double array[MAX+2][MAX+2]){
    int n;
    int k=0;
    int j=0;
    double vortex,vortexsum,mag,magsum;//无效的变量
    double restoremag,restoreavg;
    struct exp imp;
    struct record store[MAX*MAX*STEP];
    struct record statistictime(double array[MAX+2][MAX+2],struct record store[MAX*MAX*STEP],struct exp imp,int times,int k);
    printf("输入循环次数（小于%d）：",5000*MAX*MAX*STEP);
    scanf("%d",&n);
    statisticcount(array,&vortex,&vortexsum,&mag,&magsum);
    store[MAX*MAX*STEP-1].avg=vortex;
    store[MAX*MAX*STEP-1].mag=mag;
    while(k<n){
        imp=cal(array);
        k++;
        store[j++]=statistictime(array,store,imp,j-1,k);
        //print(array,vortex);
        if(k%(MAX*MAX*STEP)==0){
            //print(array,vortex);
            jilu[k/(MAX*MAX*STEP)].mag=0.0;
            jilu[k/(MAX*MAX*STEP)].avg=0.0;
            for(j=0;j<MAX*MAX*STEP;++j){
                restoremag+=store[j].mag;
                restoreavg+=store[j].avg;
                if((j+1)%(MAX*STEP)==0){
                    jilu[k/(MAX*MAX*STEP)].mag+=restoremag/(MAX*STEP);
                    jilu[k/(MAX*MAX*STEP)].avg+=restoreavg/(MAX*STEP);
                    restoremag=0.0;
                    restoreavg=0.0;
                }
            }
            jilu[k/(MAX*MAX*STEP)].mag=jilu[k/(MAX*MAX*STEP)].mag/MAX;
            jilu[k/(MAX*MAX*STEP)].avg=jilu[k/(MAX*MAX*STEP)].avg/MAX;
            jilu[k/(MAX*MAX*STEP)].number=k/(MAX*MAX*STEP)+1;
            j=0;
        }
    }
    //shellsort(jilu,n/(MAX*MAX*STEP));		//对时间演化的值排序
    printf(" 步数               涡旋平均值            磁化平均值            大小顺序\n");
    for(k=0;k<(n/(MAX*MAX*STEP));k++){		//打印随时间演化的值
        static double i=0.0,j=0.0;
        i+=jilu[k].avg;
        j+=jilu[k].mag;
        printf("   %.1f          %.16f         %.16f      %3d\n",(double)(jilu[k].number-1)*STEP,i/(k+1),j/(k+1),k+1);
    }
}


//重复多次求平均值//
void counttrace(double array[MAX+2][MAX+2]){
    int s,n;
    int q;
    double vortex,vortexsum=0.0;
    double mag,magsum=0.0;
    printf("输入循环次数：");
    scanf("%d",&n);
    printf("输入重复次数（小于1000次）：");
    scanf("%d",&s);
    for(q=0;q<s;++q){
        int k=0;
        while(k<n){
            cal(array);
            k++;
        }
        statisticcount(array,&vortex,&vortexsum,&mag,&magsum);
        print(array,vortex);
    }
    shellsort(jilu,s);			//对重复演化的值排序
    printf(" 出现顺序               涡旋平均值            磁化平均值            大小顺序\n");
    for(q=0;q<s;q++)			//打印重复演化的值
        printf("   %3d          %.16f         %.16f      %3d\n",jilu[q].number,jilu[q].avg,jilu[q].mag,q+1);
    printf("\n");
    printf("涡旋总平均：%.16f\n",vortexsum/s);
    printf("磁化总平均：%.16f\n",magsum/s);
}


//空间关联函数的数值模拟//
void relation(double array[MAX+2][MAX+2]){
    int n,j;
    int k=0,i=0;
    int x,y,z=0;
    int times[100];
    struct exp imp;
    struct val prin[100];
    double relationsum[100];
    double sum[100];
    double resum[100];
    double m[100];
    double relationcount(double array[MAX+2][MAX+2],struct val prin,int *times);
    double relationtime(double array[MAX+2][MAX+2],struct exp imp,struct val prin);
    void shellsortdouble(struct val s[],int len);
    printf("输入循环次数(每%d次输出一个值)：",MAX*MAX*STEP);scanf("%d",&n);
    for(x=0;x<TIME+1;++x){
        for(y=0;y<TIME+1;++y){
            if(sqrt(x*x+y*y)<=TIME&&(x+y)!=0){
                prin[z].x=x;
                prin[z].y=y;
                prin[z++].value=sqrt(x*x+y*y);
            }
        }
    }
    shellsortdouble(prin,z);
    for(j=0;j<z;++j){
        times[j]=0;
        sum[j]=0.0;
        m[j]=0.0;
        relationsum[j]=relationcount(array,prin[j],&times[j]);
    }
    printf("距离 ");
    for(j=0;j<z;++j){
        if(j==z-1||prin[j].value!=prin[j+1].value)
            printf("%.4f  ",prin[j].value);
    }
    printf("\n\n");
    while(k<n){
        imp=cal(array);
        for(j=0;j<z;++j){
            relationsum[j]+=relationtime(array,imp,prin[j]);
            resum[j]+=relationsum[j];
            /*double o;
             int a;
             o=relationcount(array,prin[j],&a);
             if(fabs(o-relationsum[j])>0.00001)
             printf("%d ",k);*/
            if((k+1)%(MAX*STEP)==0){
                sum[j]+=resum[j]/(MAX*STEP);
                resum[j]=0;
            }
        }
        k++;
        if(k%(MAX*MAX*STEP)==0){
            ++i;
            int time[100];
            double strm[100];
            for(j=0;j<z;++j){
                m[j]+=sum[j]/MAX;
                time[j]=times[j];
                strm[j]=m[j];
            }
            printf(" %d   ",i);
            for(j=0;j<z;++j){
                if(j==z-1||prin[j].value!=prin[j+1].value)
                    printf("%.4f  ",strm[j]/(i*time[j]));
                else{
                    strm[j+1]+=strm[j];
                    time[j+1]+=time[j];
                }
            }
            printf("\n");
            for(j=0;j<z;++j){
                sum[j]=0.0;
            }
        }
    }
}


//排序//
void shellsortdouble(struct val s[],int len){
    int i,j,gap;
    struct val temp;
    for(gap=len/2;gap>0;gap/=2)
        for(i=gap;i<len;i+=1){
            for(j=i-gap;j>=0&&s[j].value>s[j+gap].value;j-=gap){
                temp=s[j+gap];
                s[j+gap]=s[j];
                s[j]=temp;
            }
        }
}


//统计全局距离为distance的空间关联函数值//
double relationcount(double array[MAX+2][MAX+2],struct val prin,int *times){
    int i,j;
    double sum=0.0;
    if(prin.x==0||prin.y==0){
        for(i=1;i<MAX+1;++i){
            for(j=1;j<MAX+1;++j){
                if(i+prin.x<MAX+1&&j+prin.y<MAX+1){
                    sum+=cos(array[i+prin.x][j+prin.y]-array[i][j]);
                    (*times)++;
                }
            }
        }
    }else{
        for(i=1;i<MAX+1;++i){
            for(j=1;j<MAX+1;++j){
                if(i+prin.x<MAX+1&&j+prin.y<MAX+1){
                    sum+=cos(array[i+prin.x][j+prin.y]-array[i][j]);
                    (*times)++;
                }
                if(i+prin.x<MAX+1&&j-prin.y>0){
                    sum+=cos(array[i+prin.x][j-prin.y]-array[i][j]);
                    (*times)++;
                }
            }
        }
    }
    return sum;
}


//统计变化一次后距离为distance的空间关联函数值改变//
double relationtime(double array[MAX+2][MAX+2],struct exp imp,struct val prin){
    double ret=0.0;
    if(prin.x==0||prin.y==0){
        if(imp.angle1!=imp.angle2){
            if(imp.x+prin.x<MAX+1&&imp.y+prin.y<MAX+1)
                ret+=cos(array[imp.x+prin.x][imp.y+prin.y]-imp.angle2)-cos(array[imp.x+prin.x][imp.y+prin.y]-imp.angle1);
            if(imp.x-prin.x>0&&imp.y-prin.y>0)
                ret+=cos(array[imp.x-prin.x][imp.y-prin.y]-imp.angle2)-cos(array[imp.x-prin.x][imp.y-prin.y]-imp.angle1);
        }
    }else{
        if(imp.angle1!=imp.angle2){
            if(imp.x+prin.x<MAX+1&&imp.y+prin.y<MAX+1)
                ret+=cos(array[imp.x+prin.x][imp.y+prin.y]-imp.angle2)-cos(array[imp.x+prin.x][imp.y+prin.y]-imp.angle1);
            if(imp.x-prin.x>0&&imp.y-prin.y>0)
                ret+=cos(array[imp.x-prin.x][imp.y-prin.y]-imp.angle2)-cos(array[imp.x-prin.x][imp.y-prin.y]-imp.angle1);
            if(imp.x+prin.x<MAX+1&&imp.y-prin.y>0)
                ret+=cos(array[imp.x+prin.x][imp.y-prin.y]-imp.angle2)-cos(array[imp.x+prin.x][imp.y-prin.y]-imp.angle1);
            if(imp.x-prin.x>0&&imp.y+prin.y<MAX+1)
                ret+=cos(array[imp.x-prin.x][imp.y+prin.y]-imp.angle2)-cos(array[imp.x-prin.x][imp.y+prin.y]-imp.angle1);
        }
    }
    return ret;
}


//计算能量变化，判断是否改变点的状态//
struct exp cal(double array[MAX+2][MAX+2]){
    const double pi=acos(-1.0);
    double H,H1;
    int x,y;
    double a,p;
    struct exp imp;
    x=rand()%MAX+1;
    y=rand()%MAX+1;
    H1=-J*(cos(array[x][y]-array[x][y-1])+
           cos(array[x][y]-array[x][y+1])+
           cos(array[x][y]-array[x+1][y])+
           cos(array[x][y]-array[x-1][y]));
    a=array[x][y];
    array[x][y]=(rand()/(RAND_MAX+1.0))*pi*2;
    H=-J*(cos(array[x][y]-array[x][y-1])+
          cos(array[x][y]-array[x][y+1])+
          cos(array[x][y]-array[x+1][y])+
          cos(array[x][y]-array[x-1][y]))-H1;
    if(H>0){
        p=rand()/(RAND_MAX+1.0);
        if(p<=exp(-H/KT)){		//在这里记录突变点的信息
            borders(array,x,y);
        }else{
            array[x][y]=a;
        }
        imp.p=p;
    }else{
        borders(array,x,y);
        imp.p=1.0;
    }
    imp.x=x;
    imp.y=y;
    imp.p1=exp(-H/KT);
    imp.angle1=a;
    imp.angle2=array[x][y];
    return imp;
}


//打印函数//
void print(double array[MAX+2][MAX+2],double vortex){
    int i,j;
    int n=0;
    for(i=1;i<MAX+1;++i){
        for(j=1;j<MAX+1;++j){
            printf("%.2f ",array[i][j]);
            if(array[i][j]>0.0)
                n++;
        }
        printf("\n");
    }
    printf("非零数量：%d\n",n);
    printf("涡旋平均值：%.16f\n",vortex);
}


//多次重复统计旋涡平均值与磁化平均值的函数//
void statisticcount(double array[MAX+2][MAX+2],double *vortex,double *vortexsum,double *mag,double *magsum){
    int i,j;
    int m;
    static int zongshu=0;
    double sum=0.0,retu=0.0;
    double a[4];
    double sz=0.0;
    const double pi=acos(-1.0);
    for(i=1;i<MAX+1;++i){
        for(j=1;j<MAX+1;++j){
            sz+=cos(array[i][j]);
            a[0]=array[i][j+1]-array[i-1][j];
            a[1]=array[i+1][j]-array[i][j+1];
            a[2]=array[i][j-1]-array[i+1][j];
            a[3]=array[i-1][j]-array[i][j-1];
            for(m=0;m<4;++m){
                if(a[m]<-pi){
                    a[m]=a[m]+2*pi;
                }else if(a[m]>pi){
                    a[m]=a[m]-2*pi;
                }
                retu=retu+a[m];
            }
            sum+=fabs(retu);
            retu=0.0;
        }
    }
    *mag=sz/(MAX*MAX);
    *magsum+=*mag;
    jilu[zongshu].mag=*mag;
    *vortex=sum/(MAX*MAX*2*pi);
    *vortexsum+=*vortex;
    jilu[zongshu].avg=*vortex;
    jilu[zongshu].number=zongshu+1;
    zongshu++;
}


//随时间统计旋涡平均值与磁化平均值的函数//
struct record statistictime(double array[MAX+2][MAX+2],struct record store[MAX*MAX*STEP],struct exp imp,int times,int k){
    int i,j;
    int m;
    double sum=0.0,sum1=0.0,retu=0.0,retu1=0.0;
    double a[4],a1[4];
    double sz=0.0,sz1=0.0;
    struct record ret;
    const double pi=acos(-1.0);
    if(times>0){
        ret.avg=store[times-1].avg;
        ret.mag=store[times-1].mag;
    }else{
        ret.avg=store[MAX*MAX*STEP-1].avg;
        ret.mag=store[MAX*MAX*STEP-1].mag;
    }
    if(imp.angle1!=imp.angle2){
        for(i=imp.x-1;i<=imp.x+1;++i){
            for(j=imp.y-1;j<=imp.y+1;++j){
                if(i>0&&i<MAX+1&&j>0&&j<MAX+1&&
                   ((i==imp.x-1&&j==imp.y)||(i==imp.x+1&&j==imp.y)||
                    (i==imp.x&&j==imp.y-1)||(i==imp.x&&j==imp.y+1))){
                       a[0]=array[i][j+1]-array[i-1][j];
                       a[1]=array[i+1][j]-array[i][j+1];
                       a[2]=array[i][j-1]-array[i+1][j];
                       a[3]=array[i-1][j]-array[i][j-1];
                       
                       array[imp.x][imp.y]=imp.angle1;
                       borders(array,imp.x,imp.y);
                       a1[0]=array[i][j+1]-array[i-1][j];
                       a1[1]=array[i+1][j]-array[i][j+1];
                       a1[2]=array[i][j-1]-array[i+1][j];
                       a1[3]=array[i-1][j]-array[i][j-1];
                       array[imp.x][imp.y]=imp.angle2;
                       borders(array,imp.x,imp.y);
                       
                       for(m=0;m<4;++m){
                           if(a[m]<-pi){
                               a[m]=a[m]+2*pi;
                           }else if(a[m]>pi){
                               a[m]=a[m]-2*pi;
                           }
                           retu=retu+a[m];
                           
                           if(a1[m]<-pi){
                               a1[m]=a1[m]+2*pi;
                           }else if(a1[m]>pi){
                               a1[m]=a1[m]-2*pi;
                           }
                           retu1=retu1+a1[m];
                       }
                       sum+=fabs(retu);
                       sum1+=fabs(retu1);
                       retu=0.0;retu1=0.0;
                   }
                if(i==imp.x&&j==imp.y){
                    sz+=cos(array[i][j]);
                    
                    array[imp.x][imp.y]=imp.angle1;
                    sz1+=cos(array[i][j]);
                    array[imp.x][imp.y]=imp.angle2;
                }
            }
        }
        ret.mag=ret.mag+sz/(MAX*MAX)-sz1/(MAX*MAX);
        ret.avg=ret.avg+sum/(MAX*MAX*2*pi)-sum1/(MAX*MAX*2*pi);
        //double o;
        //o=check(array);
        //if(fabs(o-ret.mag)>0.001){
        //	printf("%d ",k);
        //}
    }
    ret.number=times;
    return ret;
}


//排序函数//
void shellsort(struct record s[],int len){
    int i,j,gap;
    struct record temp;
    for(gap=len/2;gap>0;gap/=2)
        for(i=gap;i<len;i+=1){
            for(j=i-gap;j>=0&&s[j].avg>s[j+gap].avg;j-=gap){
                temp=s[j+gap];
                s[j+gap]=s[j];
                s[j]=temp;
            }
        }
}


//周期性边界条件//
void borders(double array[MAX+2][MAX+2],int x,int y){
    if(x==MAX)
        array[0][y]=array[MAX][y];
    if(x==1)
        array[MAX+1][y]=array[1][y];
    if(y==MAX)
        array[x][0]=array[x][MAX];
    if(y==1)
        array[x][MAX+1]=array[x][1];
}


//检查函数//
double check(double array[MAX+2][MAX+2]){
    int i,j,m;
    double a[4];
    const double pi=acos(-1.0);
    double retu=0.0,sum=0.0,sz=0.0;
    for(i=1;i<MAX+1;++i){
        for(j=1;j<MAX+1;++j){
            sz+=cos(array[i][j]);
            a[0]=array[i][j+1]-array[i-1][j];
            a[1]=array[i+1][j]-array[i][j+1];
            a[2]=array[i][j-1]-array[i+1][j];
            a[3]=array[i-1][j]-array[i][j-1];
            for(m=0;m<4;++m){
                if(a[m]<-pi){
                    a[m]=a[m]+2*pi;
                }else if(a[m]>pi){
                    a[m]=a[m]-2*pi;
                }
                retu=retu+a[m];
            }
            //if(retu>0.0001||retu<-0.0001)
            //	printf("(%d,%d) ",i,j);
            sum+=fabs(retu);
            retu=0.0;
        }
    }
    //printf("\n%.8f\n",sum/(MAX*MAX*2*pi));
    //return sum/(MAX*MAX*2*pi);
    return sz/(MAX*MAX);
}
