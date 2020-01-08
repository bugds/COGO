def makeOptionList(text):
    refseqList = text.split('\n')
    options = ''
    template = '''
                <option value="{}">{}</option>'''
    for refseq in refseqList:
        options += template.format(refseq, refseq)
    return options