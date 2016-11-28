Subroutine RMS(F,F0,Ip,Jp,rmsf)
implicit none
integer i,j,Ip,Jp
real(8) rmsf,F(Ip,Jp),F0(Ip,Jp),res(Ip,Jp)
   DO j=1,Jp
     DO i=1,Ip
     if(abs(F0(i,j))<=1e-15) then
     res(i,j)=abs(F(i,j)-F0(i,j))/(abs(F0(i,j))+1)
     else
     res(i,j)=abs(F(i,j)-F0(i,j))/abs(F0(i,j))
     end if
     end DO
    end DO
rmsf=sum(res)/(Ip*Jp)
end Subroutine RMS
