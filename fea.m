%Ben Thomson%
%A VERY limited FEA solver. Approximates the deformation and behavior of a
%beam using Hookes law

clc

%Fundamental assumptions include static analysis, linear elastic response,
%fixed end condition and point load at one end

%first, we will try a 1D approach like the example done in class

%since we are trying to model a material, k will be constant throughout
%spring constant



load_prompt = 'Please input load in N: ';
load = input(load_prompt);
length_prompt = 'Please input length of wire in m: ';
length=input(length_prompt)


%this will model the stretch of a 10mm diameter, 1m long steel wire
%total number of nodes
%num of springs = count - 1
count = 5;

A = pi*(.05^2);
E = 200*10^9;
L=length/count;
k = (E*A)/L;

%fixed end condition
f1 = 0;



%represents the spring constant for each node
node = k*[1 -1
          -1 1];
      

%the global stiffness matrix for the model
gsm = zeros(count,count);


for n=1:count-1
  
    for row=n:2:n+1
        for col=n:2:n+1
           
            if row==col
                tmp = gsm(row:row+1,col:col+1);
                if row>1
                     gsm(row:row+1,col:col+1) = node+tmp;
                else
                     gsm(row:row+1,col:col+1) = node;
                end
            
            end
          
        end
    end
end

gsm

%as of now, gsm is singular and cannnot  be solved
%we constrain this with the boundary condiiton that there is no
%displacement at the wall (u1=0)
syms f1 u2 u3 u4 u5;
arr = [f1;u2;u3;u4;u5];
u = sym('u',[count 1]);
u(1) = 0;


bcs = [f1;0;0;0;load;];

y_coords = zeros(count+1) + 1;

gsm=gsm*u;

eqs = {};
for i=1:count
    eqs{1,i}=gsm(i)==bcs(i);
end
celldisp(eqs);


[A,B] = equationsToMatrix([eqs{1,1},eqs{1,2},eqs{1,3},eqs{1,4},eqs{1,5}], [f1, u2, u3,u4,u5]);
result = linsolve(A,B)

rxn = result(1);
result(1) = 0;

result = result+length;

%"original" length

length(y_coords);
exes = 0:L/count:L;
%plot(exes,y_coords,'x')
ylim([-1,2]);
%deformed length

plot(result, zeros(count),'o');
