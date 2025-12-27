/*
reaction = Reaction('R_3OAS140')
reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
reaction.subsystem = 'Cell Envelope Biosynthesis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

ACP_c = Metabolite(
    'ACP_c',
    formula='C11H21N2O7PRS',
    name='acyl-carrier-protein',
    compartment='c')
omrsACP_c = Metabolite(
    'M3omrsACP_c',
    formula='C25H45N2O9PRS',
    name='3-Oxotetradecanoyl-acyl-carrier-protein',
    compartment='c')

reaction.add_metabolites({
    malACP_c: -1.0,
    h_c: -1.0,
    ddcaACP_c: -1.0,
    co2_c: 1.0,
    ACP_c: 1.0,
    omrsACP_c: 1.0
})

Reaction creation values:
- identifier = Comes from KEGG Reaction ID
- name = Comes from either gbk parsing or from KEGG
- subsystem = Comes from KEGG pathways
- lower and upper bound = Comes from KEGG reaction information, techinaclly </<=>/> comes and we turn them into (-1000,0)/(-1000,1000)/(0,1000)

Metobolite creation values:
- name(informal) = Comes from KEGG reaction information from this lines below (which can be found with https://rest.kegg.jp/get/rn:R01773)
DEFINITION  L-Homoserine + NAD+ <=> L-Aspartate 4-semialdehyde + NADH + H+
EQUATION    C00263 + C00003 <=> C00441 + C00004 + C00080
- formula = Comes from KEGG get request of the compound ID
- name(formal) = KEGG reaction information
- compartment = IDK, GO terms might be able to help with this since it has CC (Cellular Component) exits. Probably not though so maybe add a parser to name it according to the reaction itself, so if the reaction has 'transport' in them add the left compounds as (e)

add_metabolite values:
- Minus / Plus values = Comes from (once again) KEGG reaction information

Run the model silently first if the objective is 0 then start the checklist process which is:
- go over the reaction and calculate the net amount of each compund mostly protons
- run the model again
- Try to find gaps with gapfill(). Although it might not be feasible since it would require a 'full' model to compare against. Maybe add this as a parameter so gap filling is possible

All that remains is parsing these files and creation of a python file that houses this model. Which will also be run in this workflow with the following parameters:
- KEGG Rxn ID --> As the model.objective
- max/min ---> As the model.direction
*/