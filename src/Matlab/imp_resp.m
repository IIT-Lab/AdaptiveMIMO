function [h11 h12 h21 h22]=imp_resp(h,R)
  p0=h';
  p1xx=p0(1,:)+i*p0(2,:);
  p1xy=p0(3,:)+i*p0(4,:);
  p1yx=p0(5,:)+i*p0(6,:);
  p1yy=p0(7,:)+i*p0(8,:);

  
  t=[1:length(p1xx)]/R;
  maxpot=0;
  for p=R/2+1:R+R/2
    pot=1.0*norm(p1xx(p:R:end))^2+.0*norm(p1yy(p:R:end))^2+...
        1.0*norm(p1yx(p:R:end))^2+.0*norm(p1xy(p:R:end))^2;
    if pot>=maxpot
      maxpot=pot;
      phase=p;
    end
  end
  h11=p1xx(phase:1:end);
  h12=p1xy(phase:1:end);
  h21=p1yx(phase:1:end);
  h22=p1yy(phase:1:end);
