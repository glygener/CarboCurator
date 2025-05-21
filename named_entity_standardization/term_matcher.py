import json
import csv
from typing import Dict, List, Tuple

class GlycanMapper:
    def __init__(self):
        self.gsd_data = []
        self.substitution_map = {}
        self.gsd_index = {}
        self.term_registry = {}

    def load_resources(self, gsd_path: str, custom_terms_path: str, custom_words_path: str):
        # Load and merge GSD with custom terms
        with open(gsd_path, 'r') as f:
            self.gsd_data = json.load(f)
        with open(custom_terms_path, 'r') as f:
            self.gsd_data += json.load(f)

        # Build GSD index and registry
        for entry in self.gsd_data:
            main_term = self.normalize(entry['term'])
            self.term_registry[main_term] = {
                'display': entry['term'],
                'glytoucan': entry.get('glytoucan_ac', ''),
                'synonyms': set(self.normalize(s) for s in entry.get('synonyms', []))
            }
            all_forms = [main_term] + list(self.term_registry[main_term]['synonyms'])
            for form in all_forms:
                self.gsd_index[form] = main_term

        # Load substitution patterns with list support
        with open(custom_words_path, 'r') as f:
            custom_words = json.load(f)
        for standard, variants in custom_words.items():
            norm_standard = self.normalize(standard)
            if isinstance(variants, str):
                self.substitution_map[self.normalize(variants)] = norm_standard
            else:
                for variant in variants:
                    self.substitution_map[self.normalize(variant)] = norm_standard

    def normalize(self, term: str) -> str:
        return term.strip().lower()

    def strict_match(self, term: str) -> Tuple[str, str]:
        norm_term = self.normalize(term)
        
        # Direct match check
        if norm_term in self.gsd_index:
            return self.gsd_index[norm_term], 'gsd_strict'
        
        # Handle plural/singular variations bidirectionally
        candidates = []
        
        # Check if adding 's' makes a match (singular -> plural)
        if not norm_term.endswith('s'):
            candidates.append(norm_term + 's')
        
        # Check if removing 's' makes a match (plural -> singular)
        candidates.append(norm_term.rstrip('s'))
        
        # Check all candidate forms
        for candidate in candidates:
            if candidate in self.gsd_index:
                return self.gsd_index[candidate], 'gsd_strict'
        
        return None, ''

    def loose_match(self, term: str) -> Tuple[str, str]:
        norm_term = self.normalize(term)
        modified = norm_term
        
        # Apply all possible substitutions
        for variant, standard in self.substitution_map.items():
            modified = modified.replace(variant, standard)
        
        # Check strict match with modified term
        main_term, _ = self.strict_match(modified)
        if main_term:
            return main_term, 'gsd_loose'
        
        return None, ''

    def process_file(self, input_path: str) -> List[Dict]:
        results = []
        
        with open(input_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                original = row['glycan_entity']
                frequency = int(row['frequency'])
                
                # Match pipeline
                main_term, match_type = self.strict_match(original)
                if not main_term:
                    main_term, match_type = self.loose_match(original)
                
                # Prepare result
                result = {
                    'original': original,
                    'frequency': frequency,
                    'mapped_term': '',
                    'glytoucan': '',
                    'match_type': match_type,
                    'confidence': 2 if match_type else 0
                }
                
                if main_term:
                    result.update({
                        'mapped_term': self.term_registry[main_term]['display'],
                        'glytoucan': self.term_registry[main_term]['glytoucan']
                    })
                    # Update frequency aggregation
                    self.term_registry[main_term]['frequency'] = \
                        self.term_registry[main_term].get('frequency', 0) + frequency
                
                results.append(result)
        
        return results

    def generate_outputs(self, mappings: List[Dict]):
        # Generate TSV
        with open('gsd_mapping.tsv', 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([
                'glycan_entity', 'mapped_term', 'frequency',
                'glytoucan_id', 'mapped_to', 'confidence'
            ])
            for row in mappings:
                writer.writerow([
                    row['original'],
                    row['mapped_term'],
                    row['frequency'],
                    row['glytoucan'],
                    row['match_type'],
                    row['confidence']
                ])

        # Generate JSON
        json_output = []
        for main_term, data in self.term_registry.items():
            if 'frequency' in data:
                json_output.append({
                    'term': data['display'],
                    'frequency': data['frequency'],
                    'synonyms': sorted([data['display']] + [
                        s for s in data['synonyms']
                        if s != self.normalize(data['display'])
                    ])
                })
        
        with open('gsd_complete_match.json', 'w') as f:
            json.dump(json_output, f, indent=2, ensure_ascii=False)

if __name__ == '__main__':
    mapper = GlycanMapper()
    mapper.load_resources(
        gsd_path='gsd_v2.7.1_cleaned.json',
        custom_terms_path='./data_input/custom_terms.json',
        custom_words_path='./data_input/custom_words.json'
    )
    mappings = mapper.process_file('./data_input/frequency_glycan_entity.tsv')
    mapper.generate_outputs(mappings)