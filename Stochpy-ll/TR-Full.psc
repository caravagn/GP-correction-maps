### Full model: telegraph process plus transcription/translation

R1:
    Protein + Gene_Active > Protein + Gene_Inactive
    kbind * Protein * Gene_Active
    
R2: 
	Gene_Inactive > Gene_Active
    kunbind * Gene_Inactive

R3:
    Gene_Active > Gene_Active + mRNA
    ktranscription * Gene_Active

R4:
    mRNA > $pool
    krnadecay * mRNA
    
R5:
    mRNA > Protein + mRNA
    ktranslation * mRNA

R6:
    Protein > $pool
    kprotdecay * Protein

#Species
mRNA = 0
Gene_Active = 0
Gene_Inactive = 1
Protein = 0

# Rates
kbind = 0.01
kunbind = 0.01

ktranslation = 1
ktranscription = 1

krnadecay = 0.01
kprotdecay = 0.1

