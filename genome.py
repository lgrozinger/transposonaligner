import blastn
import records

def align(insertions, chromosomes):
    blastn.align(insertions.values(), vectors.values())
    for insertion in insertions.values():
        features = []
        for vector in vectors.values():
            features += records.get_hsp_features_from(insertion, vector)

        if features:
            best = sorted(features, key=lambda f: blastn.hsp_score(f.hsp))[-1]
            new_features = [f for f in insertion.features if f not in features]
            new_features.append(best)
            insertion.features = new_features

        insertion.transposon = best.hsp.target if features else None
        
    return None
