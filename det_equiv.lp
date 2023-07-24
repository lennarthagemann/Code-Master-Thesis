minimize
obj: 0.656 y_1_#2 - 0.1 deltastart_Flasjo_,1_#2 - 10 zup_1_#2 - 10 zdown_1_#2
subject to
c1_1: 1 xchange_1,1_#1 - 1 xchange_2,1_#1 <= 0
c2_1: -1 lreal_Flasjon___in#1 + 1 Qnomchange_Flasjon__#1 <= 0
c3_1: 1 xchange_1,1_#2 - 1 xchange_2,1_#2 <= 0
c4_1: -0.31 Qreal_Flasjon_,1_#2 + 1 w_Flasjo_,1_#2 <= 0
c5_1: -1 w_Flasjo_,1_#2 <= 0
c6_1: 1 w_Flasjo_,1_#2 - 0.0899 u_Flasjo_,1_#2 <= 0
c1: -1 lreal_Flasjon___in#1 + 1 lreal_Flasjon___out#1 = 0
c2: -1 Qnomchange_Flasjon__#1 + 1 Qnom_Flasjon___out#1 = 0
c3: -1 xchange_1,1_#1 + 1 x_1,1__out#1 = 0
c4: -1 xchange_2,1_#1 + 1 x_2,1__out#1 = 0
c5: -1 Qinflow_Flasjon__#2 + 1 lreal_Flasjon___out#2 - 1 lreal_Flasjon___in#2 + 1 Qnom_Flasjon___in#2 = 0
c6: -1 Qnomchange_Flasjon__#2 + 1 Qnom_Flasjon___out#2 = 0
c7: -1 xchange_1,1_#2 + 1 x_1,1__out#2 = 0
c8: -1 xchange_2,1_#2 + 1 x_2,1__out#2 = 0
c9: 1 Qreal_Flasjon_,1_#2 - 1 Qnom_Flasjon___in#2 = 0
c10: 1 y_1_#2 - 0.344 x_1,1__in#2 - 0.656 x_2,1__in#2 = 0
c11: 1 x_2,1__out#1 - 1 x_2,1__in#2 = 0
c12: 1 lreal_Flasjon___out#1 - 1 lreal_Flasjon___in#2 = 0
c13: 1 x_1,1__out#1 - 1 x_1,1__in#2 = 0
c14: 1 Qnom_Flasjon___out#1 - 1 Qnom_Flasjon___in#2 = 0
Bounds
lreal_Flasjon___out#1 <= 40000
xchange_1,1_#1 <= 0.0899
xchange_2,1_#1 <= 0.0899
y_1_#1 <= 0.0899
lreal_Flasjon___out#2 <= 40000
xchange_1,1_#2 <= 0.0899
xchange_2,1_#2 <= 0.0899
y_1_#2 <= 0.0899
__20_#2 <= 0
Qnom_Flasjon___out#1 <= 10
x_1,1__out#1 <= 0.0899
x_2,1__out#1 <= 0.0899
w_Flasjo_,1_#1 <= 0.0899
__20_#1 <= 100000
Qnom_Flasjon___out#2 <= 10
x_1,1__out#2 <= 0.0899
x_2,1__out#2 <= 0.0899
w_Flasjo_,1_#2 <= 0.0899
lreal_Flasjon___out#1 >= 0
zup_1_#1 >= 0
zdown_1_#1 >= 0
Qnomchange_Flasjon__#1 >= 0
xchange_1,1_#1 >= 0
xchange_2,1_#1 >= 0
Qreal_Flasjon_,1_#1 >= 0
y_1_#1 >= 0
lreal_Flasjon___out#2 >= 0
zup_1_#2 >= 0
zdown_1_#2 >= 0
Qnomchange_Flasjon__#2 >= 0
xchange_1,1_#2 >= 0
xchange_2,1_#2 >= 0
Qreal_Flasjon_,1_#2 >= 0
y_1_#2 >= 0
__20_#2 >= 0
lreal_Flasjon___in#1 = 35000
Qnom_Flasjon___in#1 = 0
x_1,1__in#1 = 0
x_2,1__in#1 = 0
Qinflow_Flasjon__#2 = 0
Qinflow_Flasjon__#1 free
lreal_Flasjon___in#2 free
Qnom_Flasjon___in#2 free
x_1,1__in#2 free
x_2,1__in#2 free
Binary
deltastart_Flasjo_,1_#1
u_Flasjo_,1_#1
deltastart_Flasjo_,1_#2
u_Flasjo_,1_#2
End
