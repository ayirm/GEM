import re


class keggParser:
    def parse_page_entries(self, entryLocation, sectionName):
        """
        Searches the entry inside the provided page, returns the part of the searched Entry \n
        entryLocation is the pageResult which is the full page of the kegg request
        """
        if not entryLocation:
            return []
        
        lines = []
        entryExists = False
        # All KEGG entries align section names. Content starts 12 spaces in.
        content_indent = "            " # 12 spaces

        for line in entryLocation.split("\n"):
            if line.startswith(sectionName):
                entryExists = True
                lines.append(line[len(sectionName):].strip())
            elif entryExists and line.startswith(content_indent):
                lines.append(line.strip())
            elif entryExists and not line.startswith(" "):
                break
        
        return lines
    
    def find_ko_numbers(self, pageResult):
        """
        ORTHOLOGY   K12524  bifunctional aspartokinase / homoserine dehydrogenase 1 [EC:2.7.2.4 1.1.1.3]
        ^ Finds this in the page and gets the K12524 part.
        """
        koLines = self.parse_page_entries(pageResult, "ORTHOLOGY")
        koID = []

        # Regex to find one or more K numbers
        koRegex = re.compile(r'(K\d{5,})')

        for line in koLines:
            matches = koRegex.findall(line)
            for ko in matches:
                koID.append(f"KO:{ko}")

        return koID

    def find_full_brite(self, pageResult):
        """
        BRITE       KEGG Orthology (KO) [BR:ecj00001]
             09100 Metabolism
              09105 Amino acid metabolism
               00260 Glycine, serine and threonine metabolism
                JW0001 (thrA)
               00270 Cysteine and methionine metabolism
                JW0001 (thrA)
               00300 Lysine biosynthesis
                JW0001 (thrA)
              09110 Biosynthesis of other secondary metabolites
               00261 Monobactam biosynthesis
                JW0001 (thrA)
            Enzymes [BR:ecj01000]
             1. Oxidoreductases
              1.1  Acting on the CH-OH group of donors
               1.1.1  With NAD+ or NADP+ as acceptor
                1.1.1.3  homoserine dehydrogenase
                 JW0001 (thrA)
             2. Transferases
              2.7  Transferring phosphorus-containing groups
               2.7.2  Phosphotransferases with a carboxy group as acceptor
                2.7.2.4  aspartate kinase
                 JW0001 (thrA)
        ^ Finds this then just returns it
        """
        briteResults = self.parse_page_entries(pageResult, "BRITE")
        if not briteResults:
            return "No BRITE section"
        
        return "\n".join(briteResults)

    def find_brite_numbers(self, pageResult):
        briteResults = self.parse_page_entries(pageResult, "BRITE")
        if not briteResults:
            return []
        
        briteID = []
        briteRegex = re.compile(r"\[(BR:[a-z0-9]+)\]")

        for line in briteResults:
            matches = briteRegex.findall(line)
            for brID in matches:
                briteID.append(brID)

        return briteID

    def find_pathways(self, pageResult):
        """
        PATHWAY     ecj00260  Glycine, serine and threonine metabolism
            ecj00261  Monobactam biosynthesis
            ecj00270  Cysteine and methionine metabolism
            ecj00300  Lysine biosynthesis
            ecj01100  Metabolic pathways
            ecj01110  Biosynthesis of secondary metabolites
            ecj01120  Microbial metabolism in diverse environments
            ecj01230  Biosynthesis of amino acids
        ^ Finds this and then seperates the ecj... parts
        
        Then returns the ecj part. No convertion to map happens here
        """
        pathwayResults = self.parse_page_entries(pageResult, "PATHWAY")
        pathwayID = []
        pathwayRegex = re.compile(r'([a-z]{3,5}\d{5})')

        for line in pathwayResults:
            # Find the first ID on the line
            match = pathwayRegex.match(line)
            if match:
                pathwayID.append(f"path:{match.group(1)}") 
        # Add 'path:' prefix for consistency, normally just map00260 works but other ones have br:, ko: etc.
        # It was purely an aesthetic choice

        return pathwayID
    
    def find_compounds(self, reactionPage):
        """
        DEFINITION  ATP + L-Aspartate <=> ADP + 4-Phospho-L-aspartate
        EQUATION    C00002 + C00049 <=> C00008 + C03082
        parts

        Then sends them back as equation(DEFINITION) and cmpd_id(EQUATION)
        """

        if not reactionPage:
            return {
                "equation" : "N/A(ERROR)",
                "cmpdID" : "N/A(ERROR)"
            }
        
        readableEQ = " ".join(self.parse_page_entries(reactionPage, "DEFINITION"))
        cmpdEQ = " ".join(self.parse_page_entries(reactionPage, "EQUATION"))

        return {
            "equation" : readableEQ if readableEQ else "N/A",
            "cmpdID" : cmpdEQ if cmpdEQ else "N/A"
        }