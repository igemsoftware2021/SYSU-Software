from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Select, PDBIO
import copy
import os
import shutil
from Bio import Align, SeqIO
import sys
from celery import Celery

app = Celery('default')
app.conf.update({
  'broker_url': 'redis://localhost:6379',
  'enable_utc': True,
  'result_backend': 'redis://localhost:6379',
  'timezone': 'Asia/Shanghai',
})

def check_dest(destFile):
    if os.path.exists(destFile):
        os.remove(destFile)

@app.task
def long_task(path, fa, id):
    print('start4')
    result_dict = {}
    result_path = os.path.join(path, 'results')

    try:
        destFile = os.path.join(path, 'modelQ.dat')
        new = path 
        result_path_dict = {}
        model_path = os.path.join(result_path, 'model')
        if not os.path.exists(destFile):
            os.system("../run_pipeline.sh %s %s"%(path+'/'+fa, result_path))
            for i in range(5):
                old = os.path.join(model_path, 'model_' + str(i+1)+'.crderr.pdb')
                destFile = os.path.join(new, 'model_' + str(i+1)+'.crderr.pdb')
                shutil.move(old, new)
                result_path_dict['model_' + str(i+1) + '_file'] = 'model_' + str(i+1)+'.crderr.pdb'
        
            old = os.path.join(model_path, 'modelQ.dat')
            destFile = os.path.join(new, 'modelQ.dat')
            shutil.move(old, new)

        newTMScoreFile = os.path.join(new, 'modelQsummary.dat')
        outScoreLines = []
        tmscore_dict = {}
        with open(destFile) as tmscore:
            for line in tmscore.readlines():
                lineContent = line.split()
                tmscore_dict[lineContent[0]] = lineContent[1]
                pathSplit = lineContent[2].split('/')
                lineContent[2] = '/'.join(pathSplit[-2:])
                lineContent[2] += '\n'
                outScoreLines.append('\t'.join(lineContent))
            tmscore.close()

        with open(newTMScoreFile,'w+') as newtmscore:
            newtmscore.writelines(outScoreLines)
            newtmscore.close()
        result_dict['TM_score'] = tmscore_dict

        result_path_dict['TM_score_file'] = 'modelQsummary.dat'
        result_dict['result_paths'] = result_path_dict
        result_dict['pdbId'] = id
    except:
        error_title = str( sys.exc_info()[0] )
        error_detail = str( sys.exc_info() )
        result_dict['error'] = error_title + ": " + error_detail

    if os.path.exists(result_path):
        os.system('rm %s -r'%result_path)
    return result_dict

def get_candidate_residual(srcPDB):
    thresh = 0.2

    # read in src structure
    pdbId = srcPDB[srcPDB.rfind('/')+1 : srcPDB.rfind('.') if srcPDB.rfind('.') != -1 else len(srcPDB)]
    srcId = 'src_'+pdbId

    parser = PDBParser(PERMISSIVE=1)
    srcStructure = parser.get_structure(srcId, srcPDB)

    # map residual to b factor
    dict = {}
    first_chain_id = srcStructure[0].get_list()[0].get_id()
    for residual in srcStructure[0][first_chain_id].get_list():
        bfactor = 0.0
        count = 0
        for atom in residual.get_list():
            bfactor += atom.get_bfactor()
            count += 1
        dict[residual.get_id()] = bfactor/count

    # choose active residual
    sorted_keys = sorted(dict.keys(), key=lambda k :dict[k])
    candidate_keys = sorted_keys[:round(len(sorted_keys) * thresh)]
    res = []
    for residualId in candidate_keys:
        res.append(srcStructure[0][first_chain_id][residualId].get_resname() + str(residualId[1]))
    return res

def getAligned(srcPDB, tarPDB):
    # alignment
    srcRecord = SeqIO.parse(srcPDB, "pdb-atom")
    tarRecord = SeqIO.parse(tarPDB, "pdb-atom")
    srcRecord = list(srcRecord)
    tarRecord = list(tarRecord)

    # aligner = Align.PairwiseAligner()
    # alignments = aligner.align(list(srcRecord)[0].seq, list(tarRecord)[0].seq)
    seq1 = str(srcRecord[0].seq)
    seq1 = ''.join(list(filter(lambda ch : not (ch == 'X'), seq1)))
    seq2 = str(tarRecord[0].seq)
    seq2 = ''.join(list(filter(lambda ch : not (ch == 'X'), seq2)))
    # alignment = alignments[0]
    # print(alignment)
    # aligned = alignment.aligned
    pos = seq2.find(seq1)

    if pos == -1:
        aligned = (-1, -1), (-1, -1)
    else:
        aligned = ((0, len(seq1)), (pos, pos+len(seq1)))
    return aligned

def alignTarToSrc(srcPDB, tarPDB):
    # alignment
    (srcAligned, tarAligned) = getAligned(srcPDB, tarPDB)
    if srcAligned[0] == -1:
        return None, None, None, None

    # create new tar structure aligned with src sturcture
    pdbId = srcPDB[srcPDB.rfind('/')+1 : srcPDB.rfind('.') if srcPDB.rfind('.') != -1 else len(srcPDB)]
    srcId = 'src_'+pdbId
    tarId = 'tar_'+pdbId

    parser = PDBParser(PERMISSIVE=1)
    srcStructure = parser.get_structure(srcId, srcPDB)
    tarStructure = parser.get_structure(tarId, tarPDB)
    newStructure = copy.deepcopy(srcStructure)
    
    new_first_chain_id = newStructure[0].get_list()[0].get_id()
    tar_first_chain_id = tarStructure[0].get_list()[0].get_id()
    newResidualList = newStructure[0][new_first_chain_id].get_list()
    tarResidualList = tarStructure[0][tar_first_chain_id].get_list()
    outResidualList = []
    # for (src_s, src_e), (tar_s, tar_e) in zip(srcAligned, tarAligned):
    #     print(src_s, src_e, tar_s, tar_e)
    for src_i, tar_i in zip(range(srcAligned[0], srcAligned[1]), range(tarAligned[0], tarAligned[1])):
        srcResidualId = newResidualList[src_i].id
        tarResidualId = tarResidualList[tar_i].id
        outResidualList.append(srcResidualId)
        for atom in newStructure[0][new_first_chain_id][srcResidualId]:
            if not tarStructure[0][tar_first_chain_id][tarResidualId].has_id(atom.get_id()):
                print(srcResidualId)
                print(tar_first_chain_id, tarResidualId, atom.get_id())
                atom.set_coord(tarStructure[0][tar_first_chain_id][tarResidualId]['CA'].get_coord())    
            else:
                atom.set_coord(tarStructure[0][tar_first_chain_id][tarResidualId][atom.get_id()].get_coord())

    # write new structure to pdb
    class ResidualSelect(Select):
        def accept_residue(self, residue):
            if residue.get_id() in outResidualList:
                return True
            else:
                return False

    newId = 'new_tar_' + pdbId
    newTarPDB = srcPDB[:srcPDB.rfind('/') + 1] + newId + '.pdb'
    io = PDBIO()
    io.set_structure(newStructure)
    io.save(newTarPDB, ResidualSelect(), preserve_atom_numbering = True)
    newTarStructure = parser.get_structure(newId, newTarPDB)
    return srcStructure, newTarStructure, newTarPDB, new_first_chain_id

def getGroupList(active_list, srcStructure, newTarStructure, distanceThresh):
    new_tar_first_chain_id = newTarStructure[0].get_list()[0].get_id()
    src_first_chain_id = srcStructure[0].get_list()[0].get_id()
    filtered_active_list = [x for x in active_list if newTarStructure[0][new_tar_first_chain_id].has_id(x)]
    filtered_list = [x.get_id()[1] for x in srcStructure[0][src_first_chain_id].get_list() if newTarStructure[0][src_first_chain_id].has_id(x.get_id())]

    group_list = []
    for residual1Id in filtered_active_list:
        residual1 = srcStructure[0][src_first_chain_id][residual1Id]
        nearResiduals = []
        for residual2 in srcStructure[0][src_first_chain_id].get_list():
            if residual2.get_id()[1] not in filtered_list or residual2.get_id()[1] is residual1Id:
                continue

            if residual1['CA'] - residual2['CA'] < distanceThresh:
                nearResiduals.append(residual2.get_id()[1])
        if len(nearResiduals) > 0:
            group_list.append((residual1Id, nearResiduals))
    return group_list

def getScore(group_list, cadResultPath, targetPath, modelPath, new_first_chain_id, number_to_active_list):
    result = {}
    for activeResidual, neighborList in group_list:
        result[number_to_active_list[activeResidual]] = {}
        try:
            os.system("bash ../cadscore_1.1662/bin/CADscore_calc.bash -D %s -t %s -m %s -i '(%s)(%s)'"%(cadResultPath, targetPath, modelPath, new_first_chain_id+str(activeResidual), ','.join([new_first_chain_id+str(x) for x in neighborList])))
            summaryPath = cadResultPath + "/targets/" + targetPath[targetPath.rfind('/')+1:] + "/models/" +  modelPath[modelPath.rfind('/')+1:] + '/summary'
            with open(summaryPath) as lines:
                for line in lines:
                    if line[0:2] == "AA":
                        result[number_to_active_list[activeResidual]]["AA"] = float(line[2:])
                    elif line[0:2] == "AS":
                        result[number_to_active_list[activeResidual]]["AS"] = float(line[2:])
            os.system('rm %s -r'%cadResultPath)
        except:
            if os.path.exists(cadResultPath):
                os.system('rm %s -r'%cadResultPath)
            return -2

    return result

def getCad(pdbPath, srcPDB, tarPDB, active_list, custom_thresh=None):
    number_to_active_list = {}
    thresh = 4.0
    if custom_thresh:
        thresh = custom_thresh
    srcPDB = os.path.join(pdbPath, srcPDB)
    tarPDB = os.path.join(pdbPath, tarPDB)
    
    srcStructure, newTarStructure, newTarPDB, new_first_chain_id = alignTarToSrc(srcPDB, tarPDB)
    
    if srcStructure is None:
        return -1

    for e in active_list:
        num = int("".join(list((filter(str.isdigit, e)))))
        number_to_active_list[num] = e
    group_list = getGroupList(number_to_active_list.keys(), srcStructure, newTarStructure, thresh)
    if len(group_list) == 0:
        return -2
    cadResultPath = os.path.join(pdbPath, 'database_for_custom')
    result = getScore(group_list, cadResultPath, srcPDB, newTarPDB, new_first_chain_id, number_to_active_list)
    return result


if __name__ == "__main__":
    print('hello')
    # print(get_candidate_residual('../example/7CQ0/7cq0.pdb'))
    # pdbPath = '../example/'
    # srcPDB = '../example/7CQ0/7cq0.pdb'
    # tarPDB = '7CQ0/results/model/model_1.crderr.pdb'
    # active_list = ['TRP124', 'ALA122', 'TRP107', 'LEU142', 'HOH222', 'VAL121', 'HOH247', 'HOH208', 'THR138', 'THR123', 'VAL137', 'GLY91', 'LEU141', 'ALA104', 'SER125', 'THR108', 'GLY106', 'THR90', 'ILE136', 'GLY126', 'ALA110', 'GLY74', 'VAL109', 'GLY145', 'GLY158', 'ASP160', 'ARG92', 'THR155', 'PHE162', 'HIS159', 'HOH211', 'SER143', 'HOH205', 'VAL105', 'GLN127', 'LEU157', 'TRP140', 'HOH220', 'HOH270', 'VAL88', 'THR147', 'GLY71', 'PHE72', 'HOH234', 'HOH207', 'HOH228', 'SER120']
    # print(getCad(pdbPath, srcPDB, tarPDB, active_list))
    print('hello')
