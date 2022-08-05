from dataclasses import dataclass

@dataclass(init=True)
class TacoData:
    """ Class for tracking taco data """
    IC: float = 0.0
    raw_data: float = 0.0
    mean: float = 0.0
    variance: float = 0.0
    start_date: float = 0.0
    end_date: float = 0.0
    fill_factor: float = 0.0
    nuNyq: float = 0.0
    numax0_flag: bool = False
    numax_var: float = 0.0
    numax_CWTMexHat: float = 0.0
    numax_Morlet: float = 0.0
    numax0: float = 0.0
    Hmax: float = 0.0
    Bmax: float = 0.0
    HBR: float = 0.0
    Pn: float = 0.0
    A1: float = 0.0
    b1: float = 0.0
    A2: float = 0.0
    b2: float = 0.0
    A3: float = 0.0
    b3: float = 0.0
    Pg: float = 0.0
    numax: float = 0.0
    sigmaEnv: float = 0.0
    lnprob: float = 0.0
    npeaks: int = 0
    DeltaNu: float = 0.0
    DeltaNu_sd: float = 0.0
    dNu02: float = 0.0
    eps_p: float = 0.0
    eps_p_sd: float = 0.0
    alpha: float = 0.0
    alpha_sd: float = 0.0
    Central_DeltaNu: float = 0.0
    Central_DeltaNu_sd: float = 0.0
    Central_eps_p: float = 0.0
    Central_eps_p_sd: float = 0.0
    Central_alpha: float = 0.0
    Central_alpha_sd: float = 0.0
    gamma0: float = 0.0
    modeIDFlag: int = 0
