"""Choose one genuine phase canary per explicitly requested source interval.

The caller's accept function independently parses and simulates the actual
whole mutant; its rejection cannot be replaced by a label or an assumption.
This helper does not prove circuits or synthesize replacement arithmetic.
"""


def generate(original,inputs,nq,groups,mode,accept):
    if mode not in ('remove','replace'):raise ValueError('Unknown correction mutation mode')
    if not 1<=len(groups)<=4:raise ValueError('Request one through four distinct erasure groups')
    labels=[];occupied_positions=set();completed=[]
    for group in groups:
        label=group['label']
        positions=list(group['positions'])
        if not isinstance(label,str) or not label or label in labels:
            raise ValueError('Missing or duplicate source correction group label')
        if not positions or len(positions)!=len(set(positions)):
            raise ValueError('Missing or duplicate correction position')
        if any(type(i) is not int or not 0<=i<len(original) or original[i][0]!='erase' for i in positions):
            raise ValueError('Correction group does not identify actual ERASE operations')
        if occupied_positions.intersection(positions):raise ValueError('Correction source groups overlap')
        labels.append(label);occupied_positions.update(positions)
        if len(positions)>40:positions=[positions[j*(len(positions)-1)//39] for j in range(40)]
        found=False
        for i in positions:
            _,target,corrections=original[i]
            for j,correction in enumerate(corrections):
                if mode=='remove':
                    variants=[(tuple(corrections[:j])+tuple(corrections[j+1:]),dict(removed=correction))]
                else:
                    variants=[]
                    occupied={target}|{q for item in corrections for q in item[1:]}
                    for q in [q for word in inputs.values() for q in word]:
                        if q in occupied:continue
                        replacement=tuple(correction[:-1])+(q,)
                        variants.append((tuple(corrections[:j])+(replacement,)+tuple(corrections[j+1:]),
                                         dict(original=correction,replacement=replacement)))
                        if len(variants)>=12:break
                for corrected,details in variants:
                    mutation=dict(details,position=i,source_segment=label)
                    changed=original[:i]+[('erase',target,corrected)]+original[i+1:]
                    if accept('erasure_correction',changed,nq,mutation,True):
                        found=True;break
                if found:break
            if found:break
        if not found:raise ValueError('No independent input-dependent phase witness for source group '+label)
        completed.append(label)
    return completed
