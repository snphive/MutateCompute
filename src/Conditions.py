from src.Str import Str
from enum import Enum


class Cond(Enum):
    # cytoplasm and nucleus had a pH of ≈7.3, mitochondria ≈8.0, ER ≈7.5 and Golgi ≈6.6
    INCELL_MAML_FX = {Str.TEMP.value: 298, Str.PH.value: 7, Str.ION_STRGTH.value: 0.05, Str.TFE.value: 0,
                      Str.STAB.value: 0, Str.CONC.value: 1}
    # cytoplasm and nucleus had a pH of ≈7.3, mitochondria ≈8.0, ER ≈7.5 and Golgi ≈6.6
    # 310.15 K == 37 degC
    INCELL_MAML_AG = {Str.TEMP.value: 310.15, Str.PH.value: 7.4, Str.ION_STRGTH.value: 0.15, Str.TFE.value: 0,
                      Str.STAB.value: 0, Str.CONC.value: 1}
    # OUTCELL_MAML = {'temp': 310.15, 'pH': 7?}
    # 298.15 K == 25 degC
    INVITRO_COND1 = {Str.TEMP.value: 298.15, Str.PH.value: 7.5, Str.ION_STRGTH.value: 0.15, Str.TFE.value: 0,
                     Str.STAB.value: 0, Str.CONC.value: 1}
    # REFERENCE FOR THE ABOVE VALUES?
    # ETC


# NOTE: The difference between the pH values of between 7.0 and 7.4 may not have any consequence on the protein
# algorithms outputs, as I assume these values are only relevant to side-chain protonation/deprotonation and the pIs
# of the amino acids are outside of this range. It may be more important when considering mitochondrial and
# golgi-resident proteins... ?
# I'm not sure how much of an impact temperature has on Agadir and on FoldX calculations. For example would it be
# important to select 298 over 310 for example...? Does it impact the stability DDG value?
#
# pKa1, pKa2 & pKa3 are dissociation constants of alpha-carboxyl group, alpha-ammonium ion, and side chain group)
# pI (isoelectric point) is pH at which a particular molecule carries no net electrical charge or is electrically
# neutral in the statistical mean.
# pKa1 for the 20 amino acids is in the range of 1.82 to 2.83. pKa2 is 8.18 to 10.60.
# AA    pI      pKa3
# A     6.0     -
# C     5.1     -
# D     2.8     3.7
# E     3.2     4.3
# F     5.5     -
# G     6.0     -
# H     7.6     6.0
# I     6.0     -
# K     9.7     10.5
# L     6.0     -
# M     5.7     -
# N     5.4     -
# P     6.3     -
# Q     5.7     -
# R     10.8    12.5
# S     5.7     -
# T     5.6     -
# V     6.0     -
# W     5.9     -
# Y     5.7     -

# Note that pH decreases with increasing temperature:
# T (degC)      pH (in pure water, no buffer, no salt)
# 0             7.47
# 25            7.00
# 50            6.63
# 100           6.14

# Concentration of Na and K ions in mammalian cells is about 150mM
# http://book.bionumbers.org/what-are-the-concentrations-of-different-ions-in-cells/

# Concentration of protein is a tricky one. Effective and actual protein concentrations? nanoMolars upto microMolars.
