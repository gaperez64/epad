class LinIneqs(LinEqs):
    def __init__(self, eqconstrs: tuple[int, ...], ineqconstrs: tuple[int, ...]):
        LinEqs.__init__(self, eqconstrs)
        # TODO: add ineq constrs too
