# GPT-2 + sentence-level context + surprisal in bits
from pathlib import Path
import logging
import pandas as pd
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = PROJECT_ROOT / "outputs" / "roamm_words_with_length_and_frequency.csv"
OUTPUT_PATH = PROJECT_ROOT / "outputs" / "roamm_words_with_surprisal.csv"

MODEL_NAME = "gpt2"  

LOG_DIR = PROJECT_ROOT / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
LOG_PATH = LOG_DIR / f"{Path(__file__).stem}.log"

logging.basicConfig(
    level = logging.INFO,
    format = "%(asctime)s - %(levelname)s - %(message)s",
    handlers = [
        logging.FileHandler(LOG_PATH, mode="w"),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)

# Load data
logger.info("Loading %s", INPUT_PATH)

words = pd.read_csv(INPUT_PATH, keep_default_na=False)
words["word_surprisal"] = float("nan")

logger.info("Loaded %d word occurrences", len(words))

# Load language model and tokenizer
logger.info("Loading model %s", MODEL_NAME)

# use_fast = True enables the fast tokenizer implementation, which is generally faster and more efficient than the standard tokenizer. It uses the Rust-based tokenizers library, which provides optimized tokenization algorithms. This can lead to significant speed improvements, especially when processing large amounts of text. Additionally, the fast tokenizer often provides better handling of special cases and edge cases in tokenization.
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME, use_fast=True)

model = AutoModelForCausalLM.from_pretrained(MODEL_NAME)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

model.to(device)
model.eval()

logger.info("Using device: %s", device)

# Compute sentence-level surprisal
for sentence_id, sentence_df in words.groupby("sentence_id", sort=False):
    sentence_df = sentence_df.sort_values("word_position")
    sentence_words = sentence_df["words"].tolist()

    # Reconstruct sentence while remembering where every ROAMM word is located in the resulting string.
    text = ""
    word_spans = []

    for word in sentence_words:
        if text:
            text += " "
        start = len(text)
        text += word
        end = len(text)

        word_spans.append((start, end))
    
    encoded = tokenizer(text, return_offsets_mapping=True, return_tensors="pt", add_special_tokens=False)
    token_ids = encoded["input_ids"].to(device)
    offsets = encoded["offset_mapping"][0].tolist()

    # GPT-2 uses its EOS token as the start context here
    bos = torch.tensor([[tokenizer.bos_token_id]], device=device)
    input_ids = torch.cat([bos, token_ids], dim=1)
    with torch.no_grad():
        logits = model(input_ids).logits

    # Position i predicts token i+1.
    log_probs = torch.log_softmax(logits[:, :-1, :], dim=-1)
    targets = input_ids[:, 1:]
    token_log_probs = log_probs.gather(2, targets.unsqueeze(-1)).squeeze(-1)[0] ## why 2? why [0]?

    # Convert natural-log probalities to surprisal in bits.
    token_surprisal_bits = (-token_log_probs / torch.log(torch.tensor(2.0, device=device))).cpu().numpy()
    word_surprisals = [0.0] * len(sentence_words)

    # Map model tokens back onto ROAMM words.
    for token_idx, (token_start, token_end) in enumerate(offsets):
        for word_idx, (word_start, word_end) in enumerate(word_spans):
            overlap = (token_end > word_start) and (token_start < word_end)
            if overlap:
                word_surprisals[word_idx] += token_surprisal_bits[token_idx]
                break
    words.loc[sentence_df.index, "word_surprisal"] = word_surprisals


## Check and save
n_missing = words["word_surprisal"].isna().sum()

if n_missing:
    logger.warning("%d words have no surprisal computed", n_missing)

logger.info(
    "Surprisal range: %.2f to %.2f bits",
    words["word_surprisal"].min(),
    words["word_surprisal"].max()
)

words["surprisal_model"] = MODEL_NAME
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
words.to_csv(OUTPUT_PATH, index=False)


logger.info("Saved %d words to %s", len(words), OUTPUT_PATH)
