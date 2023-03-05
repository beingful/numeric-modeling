class Parameters:
    __N, __M = 5, 20
    __x0 = __y0 = __t0 = 0
    __l1, __l2, __T = 5, 5, 1
    __hx, __hy, __tau = __l1/__N, __l2/__N, __T/__M

    @property
    def N(self):
        return self.__N

    @property
    def M(self):
        return self.__M

    @property
    def x0(self):
        return self.__x0

    @property
    def y0(self):
        return self.__y0

    @property
    def t0(self):
        return self.__t0

    @property
    def l1(self):
        return self.__l1

    @property
    def l2(self):
        return self.__l2

    @property
    def T(self):
        return self.__T

    @property
    def hx(self):
        return self.__hx

    @property
    def hy(self):
        return self.__hy

    @property
    def tau(self):
        return self.__tau