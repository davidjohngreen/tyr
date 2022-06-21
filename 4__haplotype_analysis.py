file = open('status-by-ID.txt', 'r')
output = open('haplotypes.txt', 'w+')
probandFile = open('proband-IDs.txt', 'r')
probands = []

for i in probandFile:
    probands.append(i.strip('\n'))


var_402 = 0
var_192 = 1
modifier = 2

ref = 0
hmz = 2


IDs = {'CAG': [], 'CAA': [], 'TAG': [], 'TAA': [], \
    'TCA': [], 'CCA': [], 'TCG': [], 'CCG': [], '': []}


def genotype_generator(n):
    if n != './.':
        genotype = [int(i) for i in n.strip().split('/')]
        return sum(genotype)

for line in file:
    if line.startswith('LPid'):
        continue
    if line.strip().split('\t')[0] in probands:
        items = line.strip().split('\t')[1:]

        """check haplotypes"""
        if genotype_generator(items[modifier]) == ref and \
        genotype_generator(items[var_192]) == hmz and \
        genotype_generator(items[var_402]) == ref:
            IDs['CAG'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == ref and \
        genotype_generator(items[var_192]) == hmz and \
        genotype_generator(items[var_402]) == hmz:
            IDs['CAA'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == hmz and \
        genotype_generator(items[var_192]) == hmz and \
        genotype_generator(items[var_402]) == ref:
            IDs['TAG'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == hmz and \
        genotype_generator(items[var_192]) == hmz and \
        genotype_generator(items[var_402]) == hmz:
            IDs['TAA'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == hmz and \
        genotype_generator(items[var_192]) == ref and \
        genotype_generator(items[var_402]) == hmz:
            IDs['TCA'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == ref and \
        genotype_generator(items[var_192]) == ref and \
        genotype_generator(items[var_402]) == hmz:
            IDs['CCA'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == hmz and \
        genotype_generator(items[var_192]) == ref and \
        genotype_generator(items[var_402]) == ref:
            IDs['TCG'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == ref and \
        genotype_generator(items[var_192]) == ref and \
        genotype_generator(items[var_402]) == ref:
            IDs['CCG'].append(line.strip().split('\t')[0])
        else:
            IDs[''].append(line.strip().split('\t')[0])




for i in IDs:
    for j in IDs[i]:
        output.write(i + '\t' + j + '\n')
