### Reduced model: telegraph process plus translation (transcription << translation)

R1:
    Protein + Gene_Active > Protein + Gene_Inactive
    kbind * Protein * Gene_Active
    
R2: 
	Gene_Inactive > Gene_Active
    kunbind * Gene_Inactive

R3:
    Gene_Active > Gene_Active + Protein
    ktranslation * Gene_Active

R4:
    Protein > $pool
    kprotdecay * Protein

#Species
Gene_Active = 1
Gene_Inactive = 0
Protein = 0

# Rates
kbind = 0.01
kunbind = 0.01

ktranslation = 10

kprotdecay = 0.1

