from pathlib import Path
import logging
import pandas as pd
from wordsprobability import get_surprisal_per_word


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_PATH = PROJECT_ROOT / "outputs" / "02_roamm_words_with_length_and_frequency.csv"
OUTPUT_PATH = PROJECT_ROOT / "outputs" / "03_roamm_words_with_surprisal.csv"

MODEL_NAME = "gpt2-small"

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
words["word_surprisal"] = float("nan")  # Initialize the column for surprisal values
logger.info("Loaded %d word occurrences", len(words))

# Compute story_level surprisal
for story_name, story_df in words.groupby("story_name", sort=False):
    logger.info("Processing story: %s (%d words, %d pages)", story_name, len(story_df), story_df["page"].nunique())

    # Preserve the original order of words in the story
    story_words = story_df["words"].astype(str).tolist()

    # Concatenate all pages into one story
    story_text = " ".join(story_words)

    # Compute corrected surprisal with GPT-2
    result = get_surprisal_per_word(text=story_text, model_name=MODEL_NAME)

    # Check word alignment
    model_words = result["word"].tolist()

    if story_words != model_words:
        mismatches = [
            (i, r, m)
            for i, (r, m) in enumerate(zip(story_words, model_words))
            if r != m
        ]
        raise ValueError(
            f"World alignment failed for story {story_name}: "
            f"ROAMM={len(story_words)},"
            f"Model={len(model_words)}"
        )

    # Assign surprisal to the corresponding word occurrences
    words.loc[story_df.index, "word_surprisal"] = result["surprisal"].astype("float64").to_numpy()
    logger.info("Finished %s: %d words aligned", story_name, len(story_words))

# check and save
n_missing = words["word_surprisal"].isna().sum()

if n_missing:
    raise ValueError(f"{n_missing} words have no surprisal computed")

logger.info(
    "Surprisal range: %.2f to %.2f bits",
    words["word_surprisal"].min(),
    words["word_surprisal"].max()
)    

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
words.to_csv(OUTPUT_PATH, index=False)
logger.info("Saved %d word occurrences with surprisal to %s", len(words), OUTPUT_PATH)