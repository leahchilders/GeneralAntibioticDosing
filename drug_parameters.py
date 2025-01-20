# 5/31/24


class parameters:
    def __init__(
        self,
        PD_model="",
        drug_name="",
        G_max=0,
        slope=0,
        kappa=0,
        p_max=0,
        p_min=0,
        MIC=0,
        hill_coeff=0,
        E_max=0,
        EC_50=0,
        half_life=0,
        p_dose=0,
        dose_period=0,
        volume=0,
        initial_bac=0,
        min_bac=0,
        tol=0,
    ):
        self.PD_model = PD_model
        self.drug_name = drug_name
        self.G_max = G_max
        self.slope = slope
        self.kappa = kappa
        self.p_max = p_max
        self.p_min = p_min
        self.MIC = MIC
        self.hill_coeff = hill_coeff
        self.E_max = E_max
        self.EC_50 = EC_50
        self.half_life = half_life
        self.p_dose = p_dose
        self.dose_period = dose_period
        self.volume = volume
        self.initial_bac = initial_bac
        self.min_bac = min_bac
        self.tol = tol
