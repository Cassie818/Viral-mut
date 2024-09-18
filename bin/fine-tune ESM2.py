import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import esm

model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t30_150M_UR50D")
model.train()

num_layers = len(model.layers)
freeze_layers = int(0.9 * num_layers)

for i in range(freeze_layers):
    for param in model.layers[i].parameters():
        param.requires_grad = False


class SpikeProteinDataset(Dataset):
    def __init__(self, sequences):
        self.sequences = sequences

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        return self.sequences[idx]


spike_sequences = [
    "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWF",
    "MRVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWF",
    "EVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWF",
    "IFLVVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWF",
]

dataset = SpikeProteinDataset(spike_sequences)
dataloader = DataLoader(dataset, batch_size=4, shuffle=True)

criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(filter(lambda p: p.requires_grad, model.parameters()), lr=1e-5)

num_epochs = 10
for epoch in range(num_epochs):
    for batch_idx, sequences in enumerate(dataloader):
        optimizer.zero_grad()

        batch_data = [("", seq) for seq in sequences]

        batch_converter = alphabet.get_batch_converter()
        _, _, batch_tokens = batch_converter(batch_data)

        outputs = model(batch_tokens)["logits"]

        target_data = [("", seq) for seq in sequences]
        _, _, target_tokens = batch_converter(target_data)

        loss = criterion(outputs.view(-1, outputs.size(-1)), target_tokens.view(-1))

        loss.backward()
        optimizer.step()

        print(f"Epoch [{epoch + 1}/{num_epochs}], Batch [{batch_idx + 1}], Loss: {loss.item():.4f}")

torch.save(model.state_dict(), "./model/esm2_spike.pth")