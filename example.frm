#include- forcer.h

Local F =
  + 1/<p1.p1>/.../<p6.p6>*Q.p3*Q.p4*vx(Q,p1,p5,p6)*vx(-p1,p2,p3)*vx(-p5,-p6,p4)*vx(-Q,-p2,-p3,-p4)
  + 1/<p1.p1>/.../<p5.p5>*vx(-Q,p2,p3)*vx(p1,-p2,p5)*vx(-p1,p4,Q)*vx(-p3,-p4,-p5)*ex(p1,p4)
;

#call Forcer(msbarexpand=4)

B ep;
Print +s;
.end
