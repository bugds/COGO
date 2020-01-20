from bottle import Bottle, request, template
import pickle

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

from fileWork import (
    pickleFile,
    unpickleFile,
    makeOptionList,
    makeRadioButtons
)

from cogCreator import (
    getSequences,
    getIsoforms,
    countGenes,
    getRefGenes,
    removeNotRef,
    blastAndCalc
)

app = Bottle(__name__)

@app.get("/")
def index():
    return template('index')

@app.post("/listbox")
def listbox():
    eMail = request.POST.get('eMail')
    analysis = request.files.analysis
    reanalysis = request.files.reanalysis

    if reanalysis and eMail:
        inputData = pickle.load(reanalysis.file)
        genesDict = countGenes(inputData['proteins'])

        data = {
            'eMail': eMail,
            'proteins': pickleFile(inputData['proteins']),
            'options': makeOptionList(genesDict),
            'blastDict': pickleFile(inputData['blastDict'])
        }

        return template('listbox', **data)

    elif analysis and eMail:
        text = analysis.file.getvalue().decode()

        proteins = getSequences(text, dict())
        proteins = getIsoforms(proteins)
        genesDict = countGenes(proteins)

        data = {
            'eMail': eMail,
            'proteins': pickleFile(proteins),
            'options': makeOptionList(genesDict),
            'blastDict': 'blastDictIsMissing'
        }

        return template('listbox', **data)

    else:
        return 'Provide your e-mail and one of the files below'

@app.post("/reference")
def reference():
    referencial = (request.forms.get('referencial')).split(';')
    eMail = request.forms.get('eMail')
    proteins = unpickleFile(request.forms.get('proteins'))
    blastDict = request.forms.get('blastDict')

    speciesGenes = getRefGenes(referencial, proteins)

    data = {
        'eMail': eMail,
        'proteins': pickleFile(proteins),
        'blastDict': blastDict,
        'referencial': pickleFile(referencial),
        'radioButtons': makeRadioButtons(speciesGenes)
    }

    return template('reference', **data)

@app.post("/calculate")
def calculate():
    eMail = request.forms.get('eMail')
    proteins = unpickleFile(request.forms.get('proteins'))
    referencial = unpickleFile(request.forms.get('referencial'))
    blastDict = request.forms.get('blastDict')

    refDict = dict()
    for species in referencial:
        refDict[species] = request.forms.get(species)

    proteins = removeNotRef(refDict, proteins)

    if blastDict == 'blastDictIsMissing':
        blastDict = None
    else:
        blastDict = unpickleFile(blastDict)
    result = blastAndCalc(referencial, eMail, proteins, blastDict)

    data = {
        'result': result
    }

    return template('result', **data)

application = app

# for local usage:
# from bottle import run
# run(app, hostname='localhost', port=9999)