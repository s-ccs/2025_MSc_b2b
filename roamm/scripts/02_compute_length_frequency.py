from pathlib import Path
import logging
import re 
import pandas as pd 
from wordfreq import zipf_frequency


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_PATH = PROJECT_ROOT / "outputs" / "roamm_words.csv"
OUTPUT_PATH = PROJECT_ROOT / "outputs" / "roamm_words_with_length_and_frequency.csv"


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


def normalize_word(word):
    word = str(word)
    # Remove punctuation only at the beginning/end.
    # Keep internal hyphens/apostrophes:
    # "(minor-planet" -> "minor-planet"
    # "don't"         -> "don't"
    word = re.sub(r"^\W+|\W+$", "", word)
    return word

logger.info("Loading %s", INPUT_PATH)
words = pd.read_csv(INPUT_PATH, keep_default_na=False)
logger.info("Loaded %d word occurrences", len(words))

words["normalized_word"] = words["words"].apply(normalize_word)
# Word length was defined as the number of alphanumeric characters after removing surrounding punctuation.
words["word_length"] = words["normalized_word"].apply(lambda word: sum(char.isalpha() for char in word))
words["frequency_zipf"] = words["normalized_word"].apply(lambda word: zipf_frequency(word, "en") if word else float("nan"))

n_missing = words["frequency_zipf"].isna().sum()
n_zero = (words["frequency_zipf"] == 0).sum()

logger.info("%d non-lexical tokens have no frequency", n_missing)
logger.info("Word length range: %d to %d", words["word_length"].min(), words["word_length"].max())
logger.info("Frequency range: %.2f to %.2f", words["frequency_zipf"].min(), words["frequency_zipf"].max())

if n_zero > 0:
    logger.warning("%d words have zero frequency", n_zero)

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
words.to_csv(OUTPUT_PATH, index=False)

logger.info("Saved output to %s", OUTPUT_PATH)
