clear
clc
n = -2: 2;

A=n .* (n >= -2 & n <= 2);
B=A;
B_rev=fliplr(B);
d=length(A)+length(B)-1;
needplust = zeros(1,d);
needplust(1:length(A))=[1:length(A)];
needplust(d-length(A)+1:d)=[length(A):-1:1];
%test=needplust test;
%length(A) test
mark=A;
needtimes=length(A)+length(B)-1;
times=1;
%Bm=1;
%A(1,Am); B(1,Bm); test


while times<=needtimes
    B_rev=circshift(B_rev, 1);
    %if times~=needtimes
       temptimes=1;
       temp=0;
    while needplust(times)>0
    temp=temp+A(1,temptimes)*B_rev(1,temptimes);
    needplust(times)=needplust(times)-1;
    temptimes=temptimes+1;
    end

    mark(times)=temp; 
    times=times+1;
end
figure(1);
stem(n,mark)

