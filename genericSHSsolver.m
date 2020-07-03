% Generic SHS Solver

%% LCFS-NP

dStateSize = 6;
cStateSize = 3;
cdSize = cStateSize*dStateSize;

syms p1 p2 p m ex real
x = sym('x',[1,cStateSize],'real');

Q = {0 p1*[x(1),x(1),0].^ex+p2*[x(1),0,0].^ex 0 0 0 0;
    [x(1)-x(2),0,0].^ex 0 p1*[x(1),x(2),x(1)-x(2)].^ex p2*[x(1),x(2),x(3)].^ex 0 0;
    0 [x(1)-x(2),x(3),0].^ex p1*[x(1),x(2),x(1)-x(2)].^ex 0 p2*[x(1),x(2),x(3)].^ex 0;
    0 [x(1)-x(2),x(3),0].^ex 0 0 0 p1*[x(1),x(2),x(1)-x(2)].^ex;
    0 0 [x(1)-x(2),0,x(3)].^ex 0 0 p1*[x(1),x(2),x(1)-x(2)].^ex;
    0 0 0 [x(1)-x(2),x(3),0].^ex p2*[x(1),x(2),x(3)].^ex p1*[x(1),x(2),x(1)-x(2)].^ex};
%% LCFS-NP-v2

dStateSize = 6;
cStateSize = 3;
cdSize = cStateSize*dStateSize;

syms p1 p2 p m ex real
x = sym('x',[1,cStateSize],'real');

Q = {0 p1*[x(1),x(1),0].^ex+p2*[x(1),0,0].^ex 0 0 0 0;
    [x(1)-x(2),0,0].^ex 0 p1*[x(1),x(2),x(1)-x(2)].^ex p2*[x(1),x(2),0].^ex 0 0;
    0 [x(1)-x(2),x(3),0].^ex p1*[x(1),x(2),x(1)-x(2)].^ex 0 p2*[x(1),x(2),x(3)].^ex 0;
    0 [x(1)-x(2),x(3),0].^ex 0 0 0 p1*[x(1),x(2),x(1)-x(2)].^ex;
    0 0 [x(1)-x(2),0,x(3)].^ex 0 0 p1*[x(1),x(2),x(1)-x(2)].^ex;
    0 0 0 [x(1)-x(2),x(3),0].^ex p2*[x(1),x(2),x(3)].^ex p1*[x(1),x(2),x(1)-x(2)].^ex};

%% LCFS-NP - Erlang-2

dStateSize = 11;
cStateSize = 3;
cdSize = cStateSize*dStateSize;

syms p1 p2 p m ex real
x = sym('x',[1,cStateSize],'real');

Q = {0 p1*[x(1),x(1),0].^ex+p2*[x(1),0,0].^ex 0 0 0 0 0 0 0 0 0;
    0 0 p1*[x(1),x(2),x(1)-x(2)].^ex p2*[x(1),x(2),x(3)].^ex 0 0 [x(1),x(2),x(3)].^ex 0 0 0 0;
    0 0 p1*[x(1),x(2),x(1)-x(2)].^ex 0 p2*[x(1),x(2),x(3)].^ex 0 0 [x(1),x(2),x(3)].^ex 0 0 0;
    0 0 0 0 0 p1*[x(1),x(2),x(1)-x(2)].^ex 0 0 [x(1),x(2),x(3)].^ex 0 0;
    0 0 0 0 0 p1*[x(1),x(2),x(1)-x(2)].^ex 0 0 0 [x(1),x(2),x(3)].^ex 0;
    0 0 0 0 p2*[x(1),x(2),x(3)].^ex p1*[x(1),x(2),x(1)-x(2)].^ex 0 0 0 0 [x(1),x(2),x(3)].^ex;
    [x(1)-x(2),0,0].^ex 0 0 0 0 0 0 p1*[x(1),x(2),x(1)-x(2)].^ex p2*[x(1),x(2),x(3)].^ex 0 0;
    0 [x(1)-x(2),x(3),0].^ex 0 0 0 0 0 p1*[x(1),x(2),x(1)-x(2)].^ex 0 p2*[x(1),x(2),x(3)].^ex 0;
    0 [x(1)-x(2),x(3),0].^ex 0 0 0 0 0 0 0 0 p1*[x(1),x(2),x(1)-x(2)].^ex;
    0 0 [x(1)-x(2),0,x(3)].^ex 0 0 0 0 0 0 0 p1*[x(1),x(2),x(1)-x(2)].^ex;
    0 0 0 [x(1)-x(2),x(3),0].^ex 0 0 0 0 0 p2*[x(1),x(2),x(3)].^ex p1*[x(1),x(2),x(1)-x(2)].^ex};

%% LCFS-P
dStateSize = 2;
cStateSize = 3;
cdSize = cStateSize*dStateSize;

syms p1 p2 p m ex real
x = sym('x',[1,cStateSize],'real');

Q = {[x(1)-x(2),x(3),0].^ex+p1*[x(1),x(1),0].^ex p2*[x(1),0,x(2)].^ex;
    p1*[x(1),x(1),0].^ex+[x(1)-x(2),x(3),0].^ex 0};
%% LCFS-P Erlang-2
dStateSize = 4;
cStateSize = 3;
cdSize = cStateSize*dStateSize;

syms p1 p2 p m ex real
x = sym('x',[1,cStateSize],'real');

Q = {p1*[x(1),x(1),0].^ex p2*[x(1),0,x(2)].^ex [x(1),x(2),x(3)].^ex 0;
    p1*[x(1),x(1),0].^ex 0 0 [x(1),x(2),x(3)].^ex;
    p1*[x(1),x(1),0].^ex+[x(1)-x(2),x(3),0].^ex p2*[x(1),0,x(2)].^ex 0 0;
    p1*[x(1),x(1),0].^ex+[x(1)-x(2),x(3),0].^ex p2*[x(1),x(2),x(3)].^ex 0 0};
%% LCFS-P Erlang-3
dStateSize = 6;
cStateSize = 3;
cdSize = cStateSize*dStateSize;

syms p1 p2 p m ex real
x = sym('x',[1,cStateSize],'real');

Q = {p1*[x(1),x(1),0].^ex p2*[x(1),0,x(2)].^ex 0 0 [x(1),x(2),x(3)].^ex 0;
    p1*[x(1),x(1),0].^ex 0 0 0 0 [x(1),x(2),x(3)].^ex;
    p1*[x(1),x(1),0].^ex+[x(1)-x(2),x(3),0].^ex p2*[x(1),0,x(2)].^ex 0 0 0 0;
    p1*[x(1),x(1),0].^ex+[x(1)-x(2),x(3),0].^ex p2*[x(1),0,x(3)].^ex 0 0 0 0;
    p1*[x(1),x(1),0].^ex p2*[x(1),0,x(2)].^ex [x(1),x(2),x(3)].^ex 0 0 0;
    p1*[x(1),x(1),0].^ex p2*[x(1),0,x(3)].^ex 0 [x(1),x(2),x(3)].^ex 0 0};
%% LCFS-NP - yates
dStateSize = 3;
cStateSize = 3;
cdSize = cStateSize*dStateSize;

syms p1 p2 p m ex real
x = sym('x',[1,cStateSize],'real');
p = p1 + p2;

Q = {0 p1*[x(1),x(1),0].^ex+p2*[x(1),x(2),0].^ex 0;
    [x(1)-x(2),x(3),0].^ex 0 p1*[x(1),x(2),x(1)-x(2)].^ex+p2*[x(1),x(2),0].^ex;
    0 [x(1)-x(2),x(3),0].^ex p1*[x(1),x(2),x(1)-x(2)].^ex+p2*[x(1),x(2),0].^ex};
%% LCFS-NP - yates - Erlang-2
dStateSize = 5;
cStateSize = 3;
cdSize = cStateSize*dStateSize;

syms p1 p2 p m ex real
x = sym('x',[1,cStateSize],'real');
p = p1 + p2;

Q = {0 p1*[x(1),x(1),0].^ex+p2*[x(1),x(2),0].^ex 0 0 0;
    0 0 p1*[x(1),x(2),x(1)-x(2)].^ex+p2*[x(1),x(2),0].^ex [x(1),x(2),x(3)].^ex 0;
    0 0 p1*[x(1),x(2),x(1)-x(2)].^ex+p2*[x(1),x(2),0].^ex 0 [x(1),x(2),x(3)].^ex;
    [x(1)-x(2),x(3),0].^ex 0 0 0 p1*[x(1),x(2),x(1)-x(2)].^ex+p2*[x(1),x(2),0].^ex;
    0 [x(1)-x(2),x(3),0].^ex 0 0 p1*[x(1),x(2),x(1)-x(2)].^ex+p2*[x(1),x(2),0].^ex};

%% LCFS-P - yates - Erlang-2
dStateSize = 2;
cStateSize = 2;
cdSize = cStateSize*dStateSize;

syms p1 p2 p m ex real
x = sym('x',[1,cStateSize],'real');

Q = {p1*[x(1),x(1)].^ex+p2*[x(1),0].^ex [x(1),x(2)].^ex;
    p1*[x(1),x(1)].^ex+p2*[x(1),0].^ex+[x(1)-x(2),0].^ex 0};
%% LCFS-P - yates - Erlang-3
dStateSize = 3;
cStateSize = 2;
cdSize = cStateSize*dStateSize;

syms p1 p2 p m ex real
x = sym('x',[1,cStateSize],'real');

Q = {p1*[x(1),x(1)].^ex+p2*[x(1),0].^ex [x(1),x(2)].^ex 0;
    p1*[x(1),x(1)].^ex+p2*[x(1),0].^ex 0 [x(1),x(2)].^ex;
    p1*[x(1),x(1)].^ex+p2*[x(1),0].^ex+[x(1)-x(2),0].^ex 0 0};

%%
q = sym([]); % Transition rate matrix without resets
B = sym(zeros(dStateSize,dStateSize,cStateSize));
A = sym(zeros(cdSize));

for i=1:dStateSize
    for k=1:dStateSize
        q(i,k) = subs(Q{i,k}(1),ex,0);
        for l=1:cStateSize
            if Q{i,k} ~= 0
                B(i,k,l) = subs(Q{i,k}(l),ex,1);
                [e,f] = coeffs(B(i,k,l),x,'All');
                for n = 1:cStateSize
                    ef = e(find(f(:)==x(n)));
                    if isempty(ef) == 0
                        A((l-1)*dStateSize + k,(n-1)*dStateSize + i) =  ef;
                    end
                end
            end
        end
    end
end
            
A = A + diag(-repmat(sum(q,2),cStateSize,1));

A00 = A(1:dStateSize,1:dStateSize);
A01 = A(1:dStateSize,dStateSize+1:end);
A10 = A(dStateSize+1:end,1:dStateSize);
A11 = A(dStateSize+1:end,dStateSize+1:end);

D = diag(sum(q,2));

R = q-D;

Rx = [[R' zeros(size(q,1),1)];ones(1,size(q,1)+1)];
sol = rref(Rx);

pis = sol(1:end-1,end);

%

C00 = A01*A11^-1*A10;

d = rref([A00-C00 -pis/m]);

solut = simplify(d(:,end));

age=simplify(sum(solut))

%% Solver
q = sym([]);
R = sym(zeros(cdSize));
b = [1 zeros(1,cStateSize-1)];
B = diag(repmat(b,1,dStateSize));

for i=1:dStateSize
    for k=1:dStateSize
        q(i,k) = subs(Q{i,k}(1),ex,0);
        if Q{i,k} ~= 0
            for l=1:cStateSize
                u = subs(Q{i,k}(l),ex,1);
                [e,f] = coeffs(u,x,'All');
                for n = 1:cStateSize
                    ef = e(find(f(:)==x(n)));
                    if isempty(ef) == 0
                        R((i-1)*cStateSize + n, (k-1)*cStateSize + l)=ef;   
end;end;end;end;end;end

d = diag(sum(q,2));
r = q-d;
rx = [[r' zeros(size(q,1),1)];ones(1,size(q,1)+1)];
sol = rref(rx);
PIs = sol(1:end-1,end);
%
PIss = repmat(PIs,1,cStateSize)';
D = diag(reshape(repmat(sum(q,2),1,cStateSize)',[],1));

rrefv = rref([(D-R)' B'*PIss(:)/m]);
v = rrefv(:,end);
V = sum(reshape(v,cStateSize,dStateSize)',1);

age=simplify(V(1))

%% Known Formulas
syms p1 p2 p m ex real

p=p1+p2; age = 1/m*( ((1+p+p^2)^2 + 2*p^3)/((1+p+p^2)*(1+p)^2) + (1+p^2/(1+p))/p1); % LCFS-W yates
p=p1+p2; age = 1/m*(1+p)/p1; % LCFS-S yates
age = 1/m * ((p1+1)/p1 + p2/(p1+1)); %hasan_LCFS-P
%hasan LCFS-NP
age = (2*p1^10*p2 + 13*p1^9*p2^2 + 19*p1^9*p2 + 2*p1^9 + 36*p1^8*p2^3 + 104*p1^8*p2^2 + 90*p1^8*p2 + 15*p1^8 + 55*p1^7*p2^4 + 236*p1^7*p2^3 + 383*p1^7*p2^2 + 253*p1^7*p2 + 48*p1^7 + 50*p1^6*p2^5 + 284*p1^6*p2^4 + 681*p1^6*p2^3 + 806*p1^6*p2^2 + 453*p1^6*p2 + 89*p1^6 + 27*p1^5*p2^6 + 191*p1^5*p2^5 + 618*p1^5*p2^4 + 1086*p1^5*p2^3 + 1066*p1^5*p2^2 + 545*p1^5*p2 + 110*p1^5 + 8*p1^4*p2^7 + 68*p1^4*p2^6 + 290*p1^4*p2^5 + 714*p1^4*p2^4 + 1056*p1^4*p2^3 + 934*p1^4*p2^2 + 460*p1^4*p2 + 98*p1^4 + p1^3*p2^8 + 10*p1^3*p2^7 + 61*p1^3*p2^6 + 219*p1^3*p2^5 + 477*p1^3*p2^4 + 650*p1^3*p2^3 + 554*p1^3*p2^2 + 277*p1^3*p2 + 64*p1^3 + 3*p1^2*p2^7 + 23*p1^2*p2^6 + 85*p1^2*p2^5 + 184*p1^2*p2^4 + 251*p1^2*p2^3 + 217*p1^2*p2^2 + 114*p1^2*p2 + 29*p1^2 + 3*p1*p2^6 + 16*p1*p2^5 + 38*p1*p2^4 + 55*p1*p2^3 + 50*p1*p2^2 + 28*p1*p2 + 8*p1 + p2^5 + 3*p2^4 + 5*p2^3 + 5*p2^2 + 3*p2 + 1)/(m*p1*(p1 + 1)*(p1^8*p2 + 6*p1^7*p2^2 + 8*p1^7*p2 + p1^7 + 15*p1^6*p2^3 + 39*p1^6*p2^2 + 32*p1^6*p2 + 6*p1^6 + 20*p1^5*p2^4 + 76*p1^5*p2^3 + 114*p1^5*p2^2 + 75*p1^5*p2 + 16*p1^5 + 15*p1^4*p2^5 + 74*p1^4*p2^4 + 162*p1^4*p2^3 + 187*p1^4*p2^2 + 110*p1^4*p2 + 25*p1^4 + 6*p1^3*p2^6 + 36*p1^3*p2^5 + 107*p1^3*p2^4 + 183*p1^3*p2^3 + 185*p1^3*p2^2 + 104*p1^3*p2 + 25*p1^3 + p1^2*p2^7 + 7*p1^2*p2^6 + 30*p1^2*p2^5 + 75*p1^2*p2^4 + 115*p1^2*p2^3 + 110*p1^2*p2^2 + 62*p1^2*p2 + 16*p1^2 + 2*p1*p2^6 + 10*p1*p2^5 + 25*p1*p2^4 + 38*p1*p2^3 + 36*p1*p2^2 + 21*p1*p2 + 6*p1 + p2^5 + 3*p2^4 + 5*p2^3 + 5*p2^2 + 3*p2 + 1));

age = (p1 + p2 + 1)^2/(m*p1); %yates_LCFS-P Erlang-2 : erlang katsayısı üst olarak beliriyor
%% LCFS-P Erlang-k (numerical)
syms p1 p2 p m ex real
p = p1 + p2;
n = 10; % erlang shape parameter

B=p2*m*(1/p*(1-(1/(1+p))^(n-1))+1);

age = 1/(m*p1)*((1+p)^n-B/(m+m*p-B));
%% LCFS-P-yates Erlang-k (numerical)
syms p1 p2 p m ex real

n = 10; % erlang shape parameter

age = (p1+p2+1)^n /(m*p1);


%% age_1 vs age_2

age2 = simplify(subs(age,[p2 m],[3-p1 1]));
a = subs(age2,p1,0.01:0.01:2.99);
plot(fliplr(a),a)

%% age vs throughput
user = 2;
k = 2; % erlang-k, for exponential use 1
age2 = simplify(subs(age,[p2 m],[(user-1)*p1 1*k]));
span = 0.1:0.1:10;
a = subs(age2,p1,span./(k*1));
plot(span,a)
%%
user = 1:9;
age2 = simplify(subs(age,[p2 m],[1-p1 2/p1]));

a = subs(age2,[p1],[1./user]);
plot(user,a)

