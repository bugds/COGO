# Make sure your data is up-to-date;
# cogCreator uses NCBI sources which are frequently updated
# Mismatches may cause errors

import sys
from io import StringIO
from Bio import SearchIO, Entrez
from Bio.Blast import NCBIWWW

rootFolder = sys.path[0]

hitlist_size = 300

class ProteinClass():
    '''Class for proteins
    '''
    def __init__(self, species, gene, refseq, good):
        '''Initialization
        :param species: Species in which proteins are synthetized
        :param gene: Coding gene
        :param refseq: Reference sequence accession number
        :param good: Boolean, if referencial protein - True
        '''
        self.species = species
        self.gene = gene
        self.refseq = refseq
        self.good = good

def getSequences(text, proteins, good=False):
    '''Gets accession numbers from text
    :param text: Text with accession numbers
    :param proteins: Dictionary for storing information about proteins
    :param good: Boolean, True for referencial proteins
    :return: Dictionary supplemented with proteins information
    '''
    if '\r\n' in text:
        splitter = '\r\n'
    else:
        splitter = '\n'

    for line in text.split(splitter):
        if (not line in proteins) and ('_' in line):
            proteins[line] = ProteinClass(None, None, line, good)
    return proteins

def getIsoforms(proteins):
    ''' Getting isoforms for all proteins in tree

    :param proteins: Dictionary for storing information about proteins
    :return: Dictionary supplemented with isoforms
    '''
    record = Entrez.read(Entrez.elink(
        db="gene",
        dbfrom="protein",
        id=','.join(proteins.keys())
    ))
    genes = [i['Id'] for i in record[0]['LinkSetDb'][0]['Link']]

    # Safe efetch usage is less than 20 uids at a time (???)
    # IDK why, but if more it sends "An error has occured"

    for i in range(0, 1 + (len(genes) // 20)):
        if i*20 > len(genes):
            efetchAndParse(genes[i*20:], proteins)
        else:
            efetchAndParse(genes[i*20:(i+1)*20], proteins)

    proteins = checkProteins(proteins)

    return proteins

def efetchAndParse(genesPart, proteins):
    record = Entrez.efetch(
        db="gene",
        rettype="gene_table",
        retmode="text",
        id=','.join(genesPart)
    )

    genesStrings = list()
    line = record.readline()
    if "An error has occured. Please try again later." in line:
        raise TypeError("Something wrong on NCBI side. Try again")
    while '[' in line:
        genesStrings.append(line)
        line = record.readline()
    genesStrings.append(line)

    i=-1
    while line:
        if 'from: ' in line:
            i+=1
            species = genesStrings[i].split('[')[1][:-2]
            gene = genesStrings[i+1].split(',')[0].split(': ')[1]
        if 'Exon table' in line:
            refseq = line.split()[-1]
            if refseq in proteins:
                proteins[refseq].gene = gene
                proteins[refseq].species = species
            elif refseq[1] == 'P':
                proteins[refseq] = ProteinClass(species, gene, refseq, False)
                # We don't think every protein, not mentioned in input file
                # is not referencial. There are 2 cases:
                # 1) there are isoforms of this protein in the input file
                #     (if gene is referencial, "goodGeneMakesGoodProtein" helps
                #      to maintain its "good"ness, other stay "bad")
                # 2) no isoforms in the input file, old ANs are deprecated
                #     ("checkProteins" helps new proteins to inherit "good"
                #      attributes from the old ones)
        line = record.readline()

    return proteins

def checkProteins(proteins):
    ''' If outdated proteins are present in the input files,
    NCBI fetches not them, but updated ones. In this situation
    "proteins" contains 2 numbers for the same protein: one old,
    with no gene and species info, but with correct "good"
    attribute, and new, with no "good" attribute, but with gene
    and species information. This function deletes old and
    supplies new with "good" attribute
    :param proteins: Dictionary for storing information about proteins
    :return: Dictionary with deleted old and supplied new proteins
    '''

    toDel = set()

    for old in proteins.values():

        if old.species == None:
            print(old.refseq)
            record = Entrez.read(Entrez.efetch(
                db="protein",
                rettype='gp',
                retmode='xml',
                id=old.refseq
            ))

            tempor = [r['GBFeature_quals'] for r in record[0]['GBSeq_feature-table'] \
                if r['GBFeature_key'] == 'CDS']
            tempor2 = [r['GBQualifier_value'] for r in next(iter(tempor)) \
                if r['GBQualifier_name'] == 'db_xref']
            gene = tempor2[0].split(':')[1]

            for new in proteins.values():
                if new.gene == gene:
                    print(old.refseq)
                    toDel.add(old.refseq)

    for old in toDel:
        del proteins[old]

    return proteins

def unifyProteins(proteins):
    for p in proteins.values():
        p.good = False
    return proteins

def countGenes(proteins):
    genesDict = dict()

    for p in proteins.values():
        genesDict[p.gene] = p.species

    speciesList = list(genesDict.values())

    return {i:speciesList.count(i) for i in speciesList}

def getRefGenes(referencial, proteins):
    speciesGenes = dict()

    for goodSpecies in referencial:
        speciesGenes[goodSpecies] = dict()
        genesSet = set(
            [p.gene for p in proteins.values() if p.species == goodSpecies]
        )
        for gene in genesSet:
            record = Entrez.efetch(
                db="gene",
                rettype='gene_table',
                retmode='text',
                id=gene
            )

            line = record.readline()

            speciesGenes[goodSpecies][gene] = line

    return speciesGenes

def removeNotRef(refDict, proteins):
    toDel = set()

    for p in proteins.values():
        if (p.species in refDict.keys()) and (p.gene != refDict[p.species]):
            toDel.add(p.refseq)

    for refseq in toDel:
        proteins.pop(refseq)

    return proteins

    return proteins

def makeGoodProteins(goodSpecies, proteins):
    for p in proteins.values():
        if p.species in goodSpecies:
            p.good = True

    return proteins

def blastSearch(query, species):
    '''Run BLAST, save results of a search to a file and return its contents
    :param query: String with accession numbers divided by paragraphs
    :param species: String with all species, against which BLAST is performed
    :param filename: Name of original fasta file for saving results of BLAST
    '''
    records = NCBIWWW.qblast(
        'blastp',
        'refseq_protein',
        query,
        entrez_query = species,
        hitlist_size = hitlist_size
    )
    return SearchIO.parse(records.getvalue(), 'blast-xml')

def createBlastDict(blast, blastDict):
    '''Create dictionary containing BLAST results
    :param blast: contents of the XML-file with BLAST results
    :return: Dictionary containing BLAST results
    '''
    for record in blast:
        if record.id not in blastDict:
            blastDict[record.id] = {}
        for hit in record:
            species = hit.description.split('[')[1][:-1]
            if not species in blastDict[record.id]:
                substrings = hit.id.split('|')
                for i in range(len(substrings)):
                    if substrings[i] == 'ref':
                        blastDict[record.id][species] = substrings[i+1]
    return blastDict

def checkBlastDict(blastDict, proteins, iteration, previous={'queries':set(),
    'species':set()}):
    '''Checks if BLAST found all species in each case
    :param filename: Name of currently explored file
    :param blastDict: Dictionary containing BLAST results
    :param proteins: Dictionary for storing information about proteins
    :param iteration: Number of additional BLAST required currently
    :returns: "blastDict" supported with BLAST results
    '''
    speciesSet = set([p.species for p in proteins.values()])
    speciesForBlast = set()
    queriesForBlast = set()
    for record in blastDict.keys():
        lostSpecies = speciesSet - set(blastDict[record].keys())
        if bool(lostSpecies):
            speciesForBlast = speciesForBlast | lostSpecies
            queriesForBlast.add(record)
    if bool(queriesForBlast):
        if (previous['queries'] == queriesForBlast) and \
           (previous['species'] == speciesForBlast):
            for s in speciesForBlast:
                for q in queriesForBlast:
                    blastDict[q][s] = 'NA'
        else:
            newBlast = blastSearch(
                '\n'.join(queriesForBlast),
                ' OR '.join(speciesForBlast)
            )
        blastDict = createBlastDict(newBlast, blastDict)
        return checkBlastDict(blastDict, proteins, iteration + 1,
        {'queries':queriesForBlast, 'species':speciesForBlast})
    return blastDict

def checkGood(blastDict, proteins):
    ''' Consists of 2 checks:
    oneGenePerSpecies - so there would be no duplication events across
    all good proteins. Prints warning
    trulyGood - checks if BLAST results confirm that all good proteins
    do find each other via reciprocal BLAST. Prints warning
    :param blastDict: Dictionary containing BLAST results
    :param proteins: Dictionary for storing information about proteins
    :return: "proteins"
    '''
    goodProteins = [p for p in proteins.values() if p.good]

    # def oneGenePerSpecies(goodProteins=goodProteins):
    speciesSet = set()
    genesSet = set()
    for p in goodProteins:
        if (p.species in speciesSet) and (not p.gene in genesSet):
            print('WARNING!!! Your good set contains paralogs in ' + p.species)
        speciesSet.add(p.species)
        genesSet.add(p.gene)

    # def trulyGood(blastDict=blastDict, proteins=proteins):
    for refseq in proteins:
        if proteins[refseq].good:
            for species in set([p.species for p in proteins.values() if p.good]):
                if blastDict[refseq][species] in proteins:
                    if not proteins[blastDict[refseq][species]].good:
                        print('WARNING!!! Good protein:')
                        print(refseq + ' does not find another good protein in ' \
                            + species)
                        print('Finds ' + blastDict[refseq][species] + ' instead of ' \
                            + str([p.refseq for p in goodProteins if p.species == species]))
                else:
                    print('WARNING!!! Good protein:')
                    print(refseq + ' does not find another good protein in ' \
                        + species)
                    print('Finds ' + blastDict[refseq][species] + ' instead of ' \
                        + str([p.refseq for p in goodProteins if p.species == species]))

def analyzeBlastDict(blastDict, proteins):
    '''Analysis of a BLAST dictionary
    :param blastDict: Dictionary containing BLAST results
    :param proteins: Dictionary for storing information about proteins
    :return: HTML-string containing analysis results
    '''
    htmlFull = ''
    gProteins = [p for p in proteins.values() if p.good]
    gSpecies = set([p.species for p in gProteins])

    for qSpecies in set([p.species for p in proteins.values()]):
        qGenes = set()
        qReverse = dict()
        qForward = dict()

        for qGene in [p.gene for p in proteins.values() if p.species == qSpecies]:
            qGenes.add(qGene)
            qReverse[qGene] = dict() # for qRefseq use qReverse.keys()
            qForward[qGene] = set()
            goodGene = False
            # Reciprocal BLAST
            for qRefseq in [p.refseq for p in proteins.values() if p.gene == qGene]:
                if not proteins[qRefseq].good:
                    qReverse[qGene][qRefseq] = set()
                    if proteins[qRefseq].species in gSpecies:
                        raise ValueError("Same species in both groups - remove one")
                    for s in gSpecies:
                        if blastDict[qRefseq][s] in [p.refseq for p in gProteins]:
                            qReverse[qGene][qRefseq].add(s)
                else:
                    goodGene = True
            # Forward BLAST
            for gRefseq in [p.refseq for p in gProteins]:
                if blastDict[gRefseq][qSpecies] in proteins:
                    if proteins[blastDict[gRefseq][qSpecies]].gene == qGene:
                        qForward[qGene].add(gRefseq)

        if not goodGene:
            htmlFull += writeHtml(
                    proteins,
                    gSpecies,
                    set([p.refseq for p in gProteins]),
                    qSpecies,
                    qGenes,
                    qReverse,
                    qForward
            )

    return htmlFull

def writeHtml(
    proteins,
    gSpecies,
    gRefseqs,
    qSpecies,
    qGenes,
    qReverse,
    qForward):
    '''Writes BLAST analysis for single gene in form of an HTML-file
    :param proteins: Dictionary for storing information about proteins
    :param gSpecies: Set of good species
    :param gRefseqs: Set of good accession numbers
    :param qSpecies: Species linked to analyzed genes
    :param qGenes: Set of analyzed genes
    :param qReverse: Dictionary of reciprocal BLAST: isoform>species
    :param qForward: Set of species found by forward BLAST
    :return: HTML-string of BLAST analysis for single species
    '''
    htmlPart = StringIO()
    htmlString = list()
    htmlString.append('<details>\n\t<summary>{}</summary>\n')
    htmlString.append('\t<details>\n\t\t<summary>&emsp;Gene id: {}</summary>\n\t\t<details>\n\t\t\t<summary>&emsp;&emsp;{} of {} referencial proteins failed forward BLAST:</summary>\n')
    htmlString.append('\t\t\t\t&emsp;&emsp;&emsp;&emsp;{} [{}]<br>\n')
    htmlString.append('\t\t</details>\n\t\t<details>\n\t\t\t<summary>&emsp;&emsp;{} of {} isoforms failed to find all referencial proteins in first hit:</summary>\n')
    htmlString.append('\t\t\t<details>\n\t\t\t\t<summary>&emsp;&emsp;&emsp;{}: {} of {} referencial species\' proteins don\'t match :</summary>\n')
    htmlString.append('\t\t\t\t\t&emsp;&emsp;&emsp;&emsp;&emsp;{}<br>\n')
    htmlString.append('\t\t\t</details>\n')
    htmlString.append('\t\t</details>\n\t</details>\n')
    htmlString.append('</details>')
    # htmlString = [line.replace(r'\n', '\n').replace(r'\t', '\t') for line in htmlString]

    htmlPart.write(htmlString[0].format(qSpecies))
    for qGene in qGenes:
        htmlPart.write(htmlString[1].format(
        qGene,
        str(len(gRefseqs) - len(qForward[qGene])),
        len(gRefseqs)
        ))
        for fail in (gRefseqs - qForward[qGene]):
            htmlPart.write(htmlString[2].format(
                fail,
                proteins[fail].species
            ))
        htmlPart.write(htmlString[3].format(
            str(len(qReverse[qGene]) \
                - len([qR for qR in qReverse[qGene].values() if qR])),
            len(qReverse[qGene])
        ))
        for isoform, success in qReverse[qGene].items():
            htmlPart.write(htmlString[4].format(
                isoform,
                str(len(gSpecies) - len(success)),
                len(gSpecies)
            ))
            for fail in (gSpecies - success):
                htmlPart.write(htmlString[5].format(
                    fail
                ))
            htmlPart.write(htmlString[6])
        htmlPart.write(htmlString[7])
    htmlPart.write(htmlString[8])
    return htmlPart.getvalue()

def blastAndCalc(referencial, eMail, proteins, blastDict):
    Entrez.email = eMail

    if not referencial:
        return "Choose at least one referencial species!"

    proteins = makeGoodProteins(referencial, proteins)
    if not blastDict:
        blast = blastSearch(
            '\n'.join([p.refseq for p in proteins.values()]),
            ' OR '.join([p.species for p in proteins.values()])
        )
        blastDict = createBlastDict(blast, dict())
        blastDict = checkBlastDict(blastDict, proteins, 0)
    checkGood(blastDict, proteins)
    return analyzeBlastDict(blastDict, proteins)



