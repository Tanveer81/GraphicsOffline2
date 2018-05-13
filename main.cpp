#include <bits/stdc++.h>
using namespace std;

class matrix{public: double a[4][4];};
string line;
ifstream myfile ("scene.txt");
matrix* identity;
matrix* globalTransform;
stack<matrix>  xformStack;
int pc[10000];
int push = 0;
int pop = 0;
readLine(double a[]){
    if (myfile.is_open()){
      getline (myfile,line);
      //cout << line << '\n';
      vector <string> tokens;
      stringstream check1(line);
      string intermediate;
      while(getline(check1, intermediate, ' '))tokens.push_back(intermediate);

      for(int i = 0; i < tokens.size(); i++){
         a[i] = atof(tokens[i].c_str());
         //cout<< a[i] << " ";
         }
    }
    else cout << "Unable to open file";
    //cout<<"\n";
}

template <size_t size_x, size_t size_y,size_t size_x1, size_t size_y1>
matrix* matrixMul(double (&a)[size_x][size_y], double (&b)[size_x1][size_y1],int r1,int c1,int r2,int c2){
    matrix* m = new matrix();
    for(int i = 0; i < 4; ++i)for(int j = 0; j < 4; ++j)m->a[i][j]=0;
    for(int i = 0; i < r1; ++i)
        for(int j = 0; j < c2; ++j){
            for(int k = 0; k < c1; ++k){
                m->a[i][j] += a[i][k] * b[k][j];
            }
        }
    return m;
}


void printMatrix(matrix* a){
    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
            cout<< " " << a->a[i][j];
            if(j == 3)
                cout << endl;
        }
    }
}

void init(matrix* m){
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            if(i==j)m->a[i][j]=1;
            else m->a[i][j]=0;
        }
    }
}

double dot(double a[], double b[]){
    double dot=0;
    for (int i=0;i<3;i++)dot+=a[i]*b[i];
    return dot;
}

void cross(double a[], double b[], double c[]){
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=-a[0]*b[2]+a[2]*b[0];///
    c[2]=a[0]*b[1]-a[1]*b[0];
}

void rotated(double x[] , double a[] , double c[] , double angle){
    double cs = cos(angle);
    double sn = sin(angle);
    double adotx= dot(a,x);
    double acrossx[3];
    cross(a,x,acrossx);
    for(int i=0;i<3;i++)c[i]=cs*x[i] + (1-cs)*adotx*a[i] + sn*acrossx[i];
}

int main () {
    ///STAGE ONE
    cout<< fixed;
    ofstream outFile("stage1.txt");outFile<< fixed;
    ofstream outFile2("stage2.txt");outFile2<< fixed;
    ofstream outFile3("stage3.txt");outFile3<< fixed;
    identity = new matrix();
    init(identity);
    xformStack.push(*identity);
    globalTransform = identity;

    double eye[3];
    double look[3];
    double up[3];
    double gluPerspective[4];
    double triangle[3][4];
    double scale[3];
    double translate[3];
    double rotatef[4];

    readLine(eye);//for(int i = 0; i < 3; i++)cout<< eye[i] <<"\n";
    readLine(look);//for(int i = 0; i < 3; i++)cout<< look[i] <<"\n";
    readLine(up);//for(int i = 0; i < 3; i++)cout<< up[i] <<"\n";
    readLine(gluPerspective);//for(int i = 0; i < 4; i++)cout<< gluPerspective[i] <<"\n";



    ///STAGE TWO MATRIX PREPARATION
// l = look - eye
    double l[3];
    for(int i=0;i<3;i++)l[i]=look[i]-eye[i];
// l.normalize()
    double nr = sqrt(l[0]*l[0]+l[1]*l[1]+l[2]*l[2]);
    for(int i=0;i<3;i++)l[i]/=nr;
// r = l X up
    double r[3];
    cross(l,up,r);
    //for(int i=0;i<3;i++)cout<<r[i]<<" ";
// r.normalize()
    nr = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    for(int i=0;i<3;i++)r[i]/=nr;
// u = r X l
    double u[3];
    cross(r,l,u);
//    nr = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
//    for(int i=0;i<3;i++)u[i]/=nr;

    matrix* viewXformT = new matrix();
    init(viewXformT);
    for(int i=0;i<3;i++)viewXformT->a[i][3] = -eye[i];
    //printMatrix(viewXformT);

    matrix* viewXformR = new matrix();
    init(viewXformR);
    for(int i=0;i<3;i++){
        viewXformR->a[0][i] = r[i];
        viewXformR->a[1][i] = u[i];
        viewXformR->a[2][i] = -l[i];
    }
    //printMatrix(viewXformR);

    ///STAGE THREE MATRIX PREPARATION
    ///(fovY),(aspectRatio),(near),(far). gluPerspective[]

    double fovY = gluPerspective[0];
    double aspectRatio = gluPerspective[1];
    double near = gluPerspective[2];
    double far = gluPerspective[3];

    double fovX = fovY*aspectRatio;
    double t = near * tan((fovY*3.14159265359)/360);
    double rq = near * tan((fovX*3.14159265359)/360);

    matrix* projection = new matrix();
    for(int i = 0; i <= 3; i++){
        for(int j = 0; j <= 3; j++){
            projection->a[i][j]=0;
        }
    }
    projection->a[0][0]=near/rq;
    projection->a[1][1]=near/t;
    projection->a[2][2]=(far+near)/(near-far);
    projection->a[2][3]= (2*far*near)/(near-far);
    projection->a[3][2]=-1;


    while ( getline (myfile,line) ){
      if (line.compare("triangle") == 0){
        //cout << line << '\n';
        for(int i = 0; i < 3; i++){
        readLine(triangle[i]);

        ///stage1
        double r[][1]={triangle[i][0],triangle[i][1],triangle[i][2],1};
        matrix* t = matrixMul(globalTransform->a,r,4,4,4,1);

        ///stage2
        double re[][1]={t->a[0][0],t->a[1][0],t->a[2][0],1};
        matrix* rt = matrixMul(viewXformR->a,viewXformT->a,4,4,4,4);
        matrix* vt = matrixMul(rt->a,re,4,4,4,1);
        //printMatrix(vt);

        ///stage3
        double re3[][1]={vt->a[0][0],vt->a[1][0],vt->a[2][0],1};
        matrix* vt3 = matrixMul(projection->a,re3,4,4,4,1);

        double tt;
        //cout<< fixed;
        for(int j=0;j<3;j++){
            ///stage1
            if(t->a[j][0] < .0000001 && t->a[j][0] > 0 )tt = 0;
            else if(t->a[j][0] > -.0000001 && t->a[j][0] < 0 )tt = 0;
            else tt=t->a[j][0];
            //cout<< setprecision(7)<<tt<<" ";
            if(j==2)outFile<<setprecision(7)<<tt;
            else outFile<<setprecision(7)<<tt<<" ";

            ///stage2
            if(vt->a[j][0] < .0000001 && vt->a[j][0] > 0 )tt = 0;
            else if(vt->a[j][0] > -.0000001 && vt->a[j][0] < 0 )tt = 0;
            else tt=vt->a[j][0];
            //cout<< setprecision(7)<<tt<<" ";
            if(j==2)outFile2<<setprecision(7)<<tt;
            else outFile2<<setprecision(7)<<tt<<" ";

            ///stage3
            if(vt3->a[j][0] < .0000001 && vt3->a[j][0] > 0 )tt = 0;
            else if(vt3->a[j][0] > -.0000001 && vt3->a[j][0] < 0 )tt = 0;
            else tt=vt3->a[j][0]/vt3->a[3][0];
            //cout<< setprecision(7)<<tt<<" ";
            if(j==2)outFile3<<setprecision(7)<<tt;
            else outFile3<<setprecision(7)<<tt<<" ";
        }//cout<<vt3->a[3][0]<<"\n";
        outFile<<"\n";
        outFile2<<"\n";
        outFile3<<"\n";
        //cout<<"\n";
        }
        outFile<<"\n";
        outFile2<<"\n";
        outFile3<<"\n";
        //cout<<"\n";
      }

      else if (line.compare("scale") == 0){
        //cout << line << '\n';
        readLine(scale);
        matrix *as = new matrix();
        init(as);
        for(int i=0;i<3;i++)as->a[i][i]=scale[i];
        globalTransform = matrixMul(globalTransform->a,as->a,4,4,4,4);
        //printMatrix(globalTransform);
        xformStack.push(*globalTransform);
        pc[push-1]++;
      }

      else if (line.compare("translate") == 0){
        //cout << line << '\n';
        readLine(translate);
        matrix *as = new matrix();
        init(as);
        for(int i=0;i<3;i++)as->a[i][3]=translate[i];
        //printMatrix(as);
        globalTransform = matrixMul(globalTransform->a,as->a,4,4,4,4);
        //printMatrix(globalTransform);
        xformStack.push(*globalTransform);
        pc[push-1]++;
      }

      else if (line.compare("rotate") == 0){
        //cout << line << '\n';
        readLine(rotatef);
        ///1.normalize a
        double ax = rotatef[1];
        double ay = rotatef[2];
        double az = rotatef[3];
        double a = sqrt(ax*ax+ay*ay+az*az);
        ax/=a;
        ay/=a;
        az/=a;
        ///2.c1=R(i,a,angle) c2=R(j,a,angle) c3=R(k,a,angle)
        double angle = (rotatef[0]*3.14159265359)/180;
        matrix *aaa = new matrix();
        matrix c = *aaa;
        init(&c);

        double c1[3],c2[3],c3[3];
        double i[3]={1,0,0};
        double j[3]={0,1,0};
        double k[3]={0,0,1};
        double A[] = {ax,ay,az};
        rotated(i,A,c1,angle);
        rotated(j,A,c2,angle);
        rotated(k,A,c3,angle);
        for(int ii=0;ii<3;ii++){
            c.a[ii][0] = c1[ii];
            c.a[ii][1] = c2[ii];
            c.a[ii][2] = c3[ii];
        }
        globalTransform = matrixMul(globalTransform->a,c.a,4,4,4,4);
        //printMatrix(globalTransform);
        xformStack.push(*globalTransform);
        pc[push-1]++;
      }

      else if (line.compare("push") == 0){
        //cout << line << '\n';
        pc[push]=0;
        push++;
        //cout<<"push"<<push;
      }

      else if (line.compare("pop") == 0){
        //cout << line << '\n';
        for(int i =0;i<pc[push-1];i++)xformStack.pop();
        globalTransform = (matrix*)&xformStack.top();
        pop++;
        //cout<<"pop"<<pc[push-1];
        push--;
      }

      else if (line.compare("end") == 0){
        //cout << line << '\n';
      }
    }
    myfile.close();
  return 0;
}




