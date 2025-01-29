import os
from transformers import T5Tokenizer, T5EncoderModel
import torch
import re

# Define device
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

# Load the tokenizer
tokenizer = T5Tokenizer.from_pretrained('/home/dm/.cache/huggingface/transformers/prot_t5_xl_half_uniref50-enc',
                                        do_lower_case=False)

# Load the model
model = T5EncoderModel.from_pretrained("/home/dm/.cache/huggingface/transformers/prot_t5_xl_half_uniref50-enc").to(device)

# Function to process sequences and generate embeddings
def process_sequence(sequence):
    sequence = " ".join(list(re.sub(r"[UZOB]", "X", sequence)))
    ids = tokenizer(sequence, add_special_tokens=True, padding="longest")
    input_ids = torch.tensor([ids['input_ids']]).to(device)
    attention_mask = torch.tensor([ids['attention_mask']]).to(device)

    with torch.no_grad():
        embedding = model(input_ids=input_ids, attention_mask=attention_mask)

    embedding = embedding.last_hidden_state.cpu().numpy()
    features = []
    for seq_num in range(len(embedding)):
        seq_len = (attention_mask[seq_num] == 1).sum()
        seq_emd = embedding[seq_num][:seq_len - 1]
        features.append(seq_emd)
        # print(features)

    return features

# Function to read sequences from a folder and save embeddings
def process_folder(folder_path, save_path):
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        base_name, ext = os.path.splitext(filename)
        if not os.path.isfile(file_path):
            continue

        with open(file_path, 'r') as file:
            lines = file.readlines()
            if len(lines) < 2:
                continue
            protein_name = lines[0].strip()[1:]
            sequence = lines[1].strip()

            try:
                features = process_sequence(sequence)

                # output_filename = f"{protein_name}.embd"
                output_filename = f"{base_name}.embd"
                output_file_path = os.path.join(save_path, output_filename)
                with open(output_file_path, 'w') as output_file:
                    for feature in features:
                        for vector in feature:
                            feature_str = ' '.join(map(str, vector.flatten()))
                            output_file.write(feature_str + '\n')
            except RuntimeError as e:
                print(f"Error processing {filename}: {e}")

            # Clear cache after processing each file to free up memory
            torch.cuda.empty_cache()

# Example usage
folder_path = '/home/dm/ExperimentPPI/Dataset/testingset/dset_60'  # Replace with your folder path
save_path = '/home/dm/ExperimentPPI/Dataset/testingset/t5'# Replace with your folder path

os.makedirs(save_path, exist_ok=True)
process_folder(folder_path, save_path)
