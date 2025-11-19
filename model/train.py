import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, Dataset
import torch.nn.functional as F
from tqdm import tqdm
from sklearn.metrics import (roc_auc_score, accuracy_score, average_precision_score,
                             precision_score, recall_score, roc_curve, precision_recall_curve)
import matplotlib.pyplot as plt

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)


class ImprovedDataset(Dataset):
    """改进的数据集类 - 正确处理TAP和Rank特征"""
    
    def __init__(self, dataframe, aaindex_path="./data/aaindex1_pca.csv"):
        self.data = dataframe.reset_index(drop=True)
        self.aaindex = pd.read_csv(aaindex_path)  # 只读取一次
        self.amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, index):
        row = self.data.iloc[index]
        pseudosequence = row['pseudosequence']
        peptide = row['Peptide']
        label = row['Label']
        tap = row['tap_prediction_score']
        rank = row['%Rank_EL']
        
        # 提取序列特征
        pseudo_onehot = self._onehot_encoding(pseudosequence, 34)
        peptide_onehot = self._onehot_encoding(peptide, 11)
        pseudo_aa = self._get_AA_features(pseudosequence, 34)
        peptide_aa = self._get_AA_features(peptide, 11)
        
        # 拼接序列特征 (修复: 两个输入都使用拼接特征)
        pseudo_features = torch.cat([pseudo_onehot, pseudo_aa], dim=1)  # [34, 42]
        peptide_features = torch.cat([peptide_onehot, peptide_aa], dim=1)  # [11, 42]
        
        # 全局特征单独处理
        global_features = torch.tensor([tap, rank], dtype=torch.float32)
        
        return (pseudo_features.float(), 
                peptide_features.float(), 
                global_features, 
                torch.tensor(label, dtype=torch.long))
    
    def _onehot_encoding(self, sequence, maxlen):
        """One-hot编码"""
        sequence = sequence.upper()[:maxlen]
        enc_seq = torch.zeros((maxlen, 20), dtype=torch.float32)
        
        for i, aa in enumerate(sequence):
            if aa in self.amino_acids:
                enc_seq[i, self.amino_acids.index(aa)] = 1
        
        return enc_seq
    
    def _get_AA_features(self, sequence, maxlen):
        """AA特征提取"""
        sequence = sequence.ljust(maxlen, 'X')[:maxlen]
        
        all_node_feats = []
        for index, aa in enumerate(sequence):
            node_feats = self.aaindex[aa].to_list()
            
            # 位置编码
            anchar = [0, maxlen]
            seq_onehot = [0, 0]
            seq_onehot[sum([index >= i for i in anchar]) - 1] = 1
            node_feats.extend(seq_onehot)
            
            all_node_feats.append(node_feats)
        
        return torch.tensor(all_node_feats, dtype=torch.float32)


class ImprovedCNN(nn.Module):
    """改进的CNN - 正确整合全局特征"""
    
    def __init__(self):
        super(ImprovedCNN, self).__init__()
        
        # 伪序列分支 (34 x 42)
        self.pseudo_branch = nn.Sequential(
            nn.Conv2d(1, 64, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(64),
            nn.ReLU(),
            nn.Conv2d(64, 32, kernel_size=2, stride=1, padding=0),
            nn.BatchNorm2d(32),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=(2, 1), stride=2)
        )
        
        # 肽段分支 (11 x 42)
        self.peptide_branch = nn.Sequential(
            nn.Conv2d(1, 32, kernel_size=2, stride=1, padding=0),
            nn.BatchNorm2d(32),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=(2, 1), stride=2)
        )
        
        # 动态计算卷积输出维度
        self.conv_output_dim = self._get_conv_output()
        
        # 全连接层 (卷积特征 + 全局特征)
        self.classifier = nn.Sequential(
            nn.Linear(self.conv_output_dim + 2, 1024),  # +2 for TAP and Rank
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(1024, 128),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(128, 2)
        )
    
    def _get_conv_output(self):
        """自动计算卷积输出维度"""
        with torch.no_grad():
            # Use correct feature dimensions: 20 (one-hot) + 22 (AA features) = 42
            x1 = torch.zeros(1, 1, 34, 42)
            x2 = torch.zeros(1, 1, 11, 42)
            x1 = self.pseudo_branch(x1)
            x2 = self.peptide_branch(x2)
            # Calculate flattened size per sample (excluding batch dimension)
            x1_flat = x1.view(1, -1).size(1)
            x2_flat = x2.view(1, -1).size(1)
            return x1_flat + x2_flat
    
    def forward(self, pseudo, peptide, global_feat):
        # 添加通道维度
        x1 = pseudo.unsqueeze(1)  # [batch, 1, 34, 42]
        x2 = peptide.unsqueeze(1)  # [batch, 1, 11, 42]
        
        # 卷积处理
        x1 = self.pseudo_branch(x1)
        x2 = self.peptide_branch(x2)
        
        # 展平
        x1 = x1.view(x1.size(0), -1)
        x2 = x2.view(x2.size(0), -1)
        
        # 拼接: 序列特征 + 全局特征
        x = torch.cat([x1, x2, global_feat], dim=1)
        
        # 分类 (不使用softmax,让CrossEntropyLoss处理)
        x = self.classifier(x)
        return x


class OriginalCNN(nn.Module):
    """保持原始架构用于对比"""
    
    def __init__(self):
        super(OriginalCNN, self).__init__()
        
        self.conv1 = nn.Sequential(
            nn.Conv2d(1, 64, kernel_size=3, stride=1, padding=0),
            nn.BatchNorm2d(64),
            nn.ReLU()
        )
        self.conv2 = nn.Sequential(
            nn.Conv2d(64, 32, kernel_size=2, stride=1, padding=0),
            nn.BatchNorm2d(32),
            nn.ReLU()
        )
        self.conv3 = nn.Sequential(
            nn.Conv2d(1, 32, kernel_size=2, stride=1, padding=0),
            nn.BatchNorm2d(32),
            nn.ReLU()
        )
        self.pool = nn.MaxPool2d(kernel_size=(2, 1), stride=2)
        
        # 动态计算维度
        self.fc_input_dim = self._get_conv_output()
        
        self.fc1 = nn.Linear(self.fc_input_dim + 2, 1024)
        self.fc2 = nn.Linear(1024, 128)
        self.fc3 = nn.Linear(128, 2)
        self.dropout = nn.Dropout(0.5)
    
    def _get_conv_output(self):
        with torch.no_grad():
            # Use correct feature dimensions: 20 (one-hot) + 22 (AA features) = 42
            x1 = torch.zeros(1, 1, 34, 42)
            x2 = torch.zeros(1, 1, 11, 42)
            x1 = self.conv1(x1)
            x1 = self.conv2(x1)
            x1 = self.pool(x1)
            x2 = self.conv3(x2)
            x2 = self.pool(x2)
            # Calculate flattened size per sample (excluding batch dimension)
            x1_flat = x1.view(1, -1).size(1)
            x2_flat = x2.view(1, -1).size(1)
            return x1_flat + x2_flat
    
    def forward(self, input1, input2, global_feat):
        x1 = input1.unsqueeze(1)
        x2 = input2.unsqueeze(1)
        
        x1 = self.conv1(x1)
        x1 = self.conv2(x1)
        x1 = self.pool(x1)
        
        x2 = self.conv3(x2)
        x2 = self.pool(x2)
        
        x1 = x1.view(x1.size(0), -1)
        x2 = x2.view(x2.size(0), -1)
        
        x = torch.cat([x1, x2, global_feat], dim=1)
        
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        x = F.relu(self.fc2(x))
        x = self.dropout(x)
        x = self.fc3(x)
        return x


def train_epoch(model, loader, criterion, optimizer, device):
    """训练一个epoch"""
    model.train()
    total_loss = 0
    all_preds, all_labels = [], []
    
    for pseudo, peptide, global_feat, labels in tqdm(loader, desc="Training", leave=False):
        pseudo = pseudo.to(device)
        peptide = peptide.to(device)
        global_feat = global_feat.to(device)
        labels = labels.to(device)
        
        optimizer.zero_grad()
        outputs = model(pseudo, peptide, global_feat)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item()
        preds = torch.argmax(outputs, dim=1)
        all_preds.extend(preds.cpu().numpy())
        all_labels.extend(labels.cpu().numpy())
    
    avg_loss = total_loss / len(loader)
    accuracy = accuracy_score(all_labels, all_preds)
    return avg_loss, accuracy


def evaluate_model(model, criterion, dataloader, device):
    """完整的评估函数"""
    model.eval()
    all_labels = []
    all_predictions = []
    all_probs = []
    running_loss = 0.0
    
    with torch.no_grad():
        for pseudo, peptide, global_feat, labels in tqdm(dataloader, desc="Evaluation", leave=False):
            pseudo = pseudo.to(device)
            peptide = peptide.to(device)
            global_feat = global_feat.to(device)
            labels = labels.to(device)
            
            outputs = model(pseudo, peptide, global_feat)
            loss = criterion(outputs, labels)
            running_loss += loss.item()
            
            # 正确计算概率
            probs = F.softmax(outputs, dim=1)[:, 1]
            predictions = torch.argmax(outputs, dim=1)
            
            all_labels.extend(labels.cpu().numpy())
            all_predictions.extend(predictions.cpu().numpy())
            all_probs.extend(probs.cpu().numpy())
    
    avg_loss = running_loss / len(dataloader)
    
    # 计算所有指标
    precision = precision_score(all_labels, all_predictions, zero_division=0)
    recall = recall_score(all_labels, all_predictions, zero_division=0)
    auc = roc_auc_score(all_labels, all_probs)
    aupr = average_precision_score(all_labels, all_probs)
    acc = accuracy_score(all_labels, all_predictions)
    
    # ROC和PR曲线数据
    fpr, tpr, _ = roc_curve(all_labels, all_probs)
    precision_curve, recall_curve, _ = precision_recall_curve(all_labels, all_probs)
    
    return {
        'loss': avg_loss,
        'auc': auc,
        'aupr': aupr,
        'acc': acc,
        'precision': precision,
        'recall': recall,
        'fpr': fpr,
        'tpr': tpr,
        'precision_curve': precision_curve,
        'recall_curve': recall_curve
    }


def plot_training_history(history):
    """绘制训练历史"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Loss
    axes[0, 0].plot(history['train_loss'], label='Train Loss')
    axes[0, 0].plot(history['val_loss'], label='Val Loss')
    axes[0, 0].set_xlabel('Epoch')
    axes[0, 0].set_ylabel('Loss')
    axes[0, 0].set_title('Training and Validation Loss')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Accuracy
    axes[0, 1].plot(history['train_acc'], label='Train Acc')
    axes[0, 1].plot(history['val_acc'], label='Val Acc')
    axes[0, 1].set_xlabel('Epoch')
    axes[0, 1].set_ylabel('Accuracy')
    axes[0, 1].set_title('Training and Validation Accuracy')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # AUC
    axes[1, 0].plot(history['val_auc'], label='Val AUC', color='green')
    axes[1, 0].set_xlabel('Epoch')
    axes[1, 0].set_ylabel('AUC')
    axes[1, 0].set_title('Validation AUC')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # AUPR
    axes[1, 1].plot(history['val_aupr'], label='Val AUPR', color='purple')
    axes[1, 1].set_xlabel('Epoch')
    axes[1, 1].set_ylabel('AUPR')
    axes[1, 1].set_title('Validation AUPR')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('training_history.png', dpi=300, bbox_inches='tight')
    plt.show()


def plot_roc_pr_curves(metrics):
    """绘制ROC和PR曲线"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # ROC Curve
    axes[0].plot(metrics['fpr'], metrics['tpr'], 
                color='darkorange', lw=2, 
                label=f"ROC curve (AUC = {metrics['auc']:.3f})")
    axes[0].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random')
    axes[0].set_xlabel('False Positive Rate', fontsize=12)
    axes[0].set_ylabel('True Positive Rate', fontsize=12)
    axes[0].set_title('Receiver Operating Characteristic', fontsize=14)
    axes[0].legend(loc='lower right')
    axes[0].grid(True, alpha=0.3)
    
    # PR Curve
    axes[1].plot(metrics['recall_curve'], metrics['precision_curve'],
                color='blue', lw=2,
                label=f"PR curve (AUPR = {metrics['aupr']:.3f})")
    axes[1].set_xlabel('Recall', fontsize=12)
    axes[1].set_ylabel('Precision', fontsize=12)
    axes[1].set_title('Precision-Recall Curve', fontsize=14)
    axes[1].legend(loc='lower left')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('roc_pr_curves.png', dpi=300, bbox_inches='tight')
    plt.show()


def train_full_pipeline(model_class, model_name, train_loader, val_loader, 
                       num_epochs=50, lr=0.0001, device=device):
    """完整的训练流程"""
    print(f"\n{'='*60}")
    print(f"Training {model_name}")
    print(f"{'='*60}\n")
    
    model = model_class().to(device)
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=lr, betas=(0.9, 0.999), weight_decay=0.001)
    
    history = {
        'train_loss': [], 'train_acc': [],
        'val_loss': [], 'val_acc': [],
        'val_auc': [], 'val_aupr': []
    }
    
    best_auc = 0
    
    for epoch in range(num_epochs):
        # 训练
        train_loss, train_acc = train_epoch(model, train_loader, criterion, optimizer, device)
        
        # 验证
        val_metrics = evaluate_model(model, criterion, val_loader, device)
        
        # 记录历史
        history['train_loss'].append(train_loss)
        history['train_acc'].append(train_acc)
        history['val_loss'].append(val_metrics['loss'])
        history['val_acc'].append(val_metrics['acc'])
        history['val_auc'].append(val_metrics['auc'])
        history['val_aupr'].append(val_metrics['aupr'])
        
        # 打印结果
        print(f"Epoch {epoch+1}/{num_epochs}")
        print(f"  Train - Loss: {train_loss:.4f}, Acc: {train_acc:.4f}")
        print(f"  Val   - Loss: {val_metrics['loss']:.4f}, Acc: {val_metrics['acc']:.4f}, "
              f"AUC: {val_metrics['auc']:.4f}, AUPR: {val_metrics['aupr']:.4f}")
        print(f"          Precision: {val_metrics['precision']:.4f}, Recall: {val_metrics['recall']:.4f}")
        
        # 保存最佳模型
        if val_metrics['auc'] > best_auc:
            best_auc = val_metrics['auc']
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'auc': best_auc,
            }, f'best_{model_name}.pth')
            print(f"  [*] Best model saved! (AUC: {best_auc:.4f})")
        print()
    
    # 最终评估
    final_metrics = evaluate_model(model, criterion, val_loader, device)
    
    return model, history, final_metrics


if __name__ == "__main__":
    # 读取数据
    print("Loading data...")
    data = pd.read_csv("data/output_binding_tap.csv")
    train_data, val_data = train_test_split(data, test_size=0.2, random_state=42)
    
    # 创建数据集
    train_dataset = ImprovedDataset(train_data)
    val_dataset = ImprovedDataset(val_data)
    
    train_loader = DataLoader(train_dataset, batch_size=256, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=256, shuffle=False)
    
    print(f"Train samples: {len(train_dataset)}")
    print(f"Val samples: {len(val_dataset)}")
    
    # 训练改进的CNN
    model, history, final_metrics = train_full_pipeline(
        ImprovedCNN, 
        "ImprovedCNN",
        train_loader, 
        val_loader,
        num_epochs=50,
        lr=0.0001
    )
    
    print("\n" + "="*60)
    print("FINAL RESULTS")
    print("="*60)
    print(f"Accuracy:  {final_metrics['acc']:.4f}")
    print(f"AUC:       {final_metrics['auc']:.4f}")
    print(f"AUPR:      {final_metrics['aupr']:.4f}")
    print(f"Precision: {final_metrics['precision']:.4f}")
    print(f"Recall:    {final_metrics['recall']:.4f}")
    
    # 绘制训练历史
    plot_training_history(history)
    
    # 绘制ROC和PR曲线
    plot_roc_pr_curves(final_metrics)
    
    print("\nTraining completed! Figures saved.")