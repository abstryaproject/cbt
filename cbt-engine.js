/*************************************************
 * 1️⃣ QUESTION COLLECTION (already prepared)
 *************************************************/

const TOTAL_QUESTIONS = 60;


/*************************************************
 * 2️⃣ SHUFFLE FUNCTION (SAFE)
 *************************************************/
function shuffleArray(array) {
  const arr = [...array];
  for (let i = arr.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [arr[i], arr[j]] = [arr[j], arr[i]];
  }
  return arr;
}


/*************************************************
 * 3️⃣ INITIALIZE EXAM
 *************************************************/
let quizData = [];
let currentQuestion = 0;
let timerInterval;
let timeLeft = examTime; // examTime must already exist
let inPreviewMode = false;
let currentAnswerPreview = 0;


/*************************************************
 * 4️⃣ START / RESTART EXAM (NEW RANDOM QUESTIONS)
 *************************************************/
function startExam() {
  quizData = shuffleArray(allQuizData)
    .slice(0, TOTAL_QUESTIONS)
    .map(q => ({
      ...q,
      options: shuffleArray(q.options),
      userAnswer: null
    }));

  currentQuestion = 0;
  timeLeft = examTime;
  inPreviewMode = false;

  quizContainer.innerHTML = '';
  resultContainer.innerHTML = '';
  resultContainer.style.display = 'block';

  prevBtn.style.display = 'inline-block';
  nextBtn.style.display = 'inline-block';
  submitBtn.style.display = 'inline-block';
  retryBtn.style.display = 'none';
  backToScoreBtn.style.display = 'none';
  homeBtn.style.display = 'none';

  timerElement.style.display = 'block';

  startTimer();
  loadQuestion();
}


/*************************************************
 * 5️⃣ TIMER
 *************************************************/
function startTimer() {
  clearInterval(timerInterval);
  timerInterval = setInterval(() => {
    timeLeft--;
    timerElement.textContent = `Time Left: ${timeLeft}s`;
    if (timeLeft <= 0) {
      clearInterval(timerInterval);
      showResults();
    }
  }, 1000);
}


/*************************************************
 * 6️⃣ LOAD QUESTION
 *************************************************/
function loadQuestion() {
  const q = quizData[currentQuestion];
  quizContainer.innerHTML = `
    <h3>Question ${currentQuestion + 1} of ${quizData.length}</h3>
    <p>${q.question}</p>
    <div class="answers">
      ${q.options.map((opt, i) => `
        <label>
          <input type="radio" name="answer" value="${i}"
            ${q.userAnswer === i ? 'checked' : ''}>
          ${opt}
        </label>
      `).join('')}
    </div>
  `;
}


/*************************************************
 * 7️⃣ SAVE ANSWER
 *************************************************/
function saveAnswer() {
  const selected = document.querySelector('input[name="answer"]:checked');
  if (selected) quizData[currentQuestion].userAnswer = Number(selected.value);
}


/*************************************************
 * 8️⃣ NAVIGATION
 *************************************************/
nextBtn.onclick = () => {
  saveAnswer();
  if (currentQuestion < quizData.length - 1) {
    currentQuestion++;
    loadQuestion();
  }
};

prevBtn.onclick = () => {
  saveAnswer();
  if (currentQuestion > 0) {
    currentQuestion--;
    loadQuestion();
  }
};

submitBtn.onclick = showResults;


/*************************************************
 * 9️⃣ SHOW RESULTS
 *************************************************/
function showResults() {
  saveAnswer();
  clearInterval(timerInterval);
  timerElement.style.display = 'none';

  let score = 0;
  quizData.forEach(q => {
    if (q.userAnswer === q.answer) score++;
  });

  const percent = (score / quizData.length) * 100;
  const grade =
    percent >= 90 ? 'A' :
    percent >= 70 ? 'B' :
    percent >= 50 ? 'C' : 'F';

  const isFail = percent < 50;

  quizContainer.innerHTML = '';
  prevBtn.style.display = 'none';
  nextBtn.style.display = 'none';
  submitBtn.style.display = 'none';

  resultContainer.innerHTML = `
    <div style="font-size:1.1rem;">
      <b>Score:</b> ${score} / ${quizData.length} (${percent.toFixed(0)}%)<br>
      <b>Grade:</b> ${grade}<br><br>

      ${!isFail ? `<button id="viewAnswersBtn">View Correct Answers & Explanation</button>` : ''}
      ${isFail ? `<button id="retryExamBtn">Retry Exam</button>` : ''}
      <button id="homeBtn">Home</button>
    </div>
  `;

  // Home button
  document.getElementById('homeBtn').onclick = () => {
    window.location.href = 'index.html';
  };

  // Retry exam (FAIL ONLY)
  if (isFail) {
    document.getElementById('retryExamBtn').onclick = startExam;
  }

  // View answers (PASS ONLY)
  if (!isFail) {
    document.getElementById('viewAnswersBtn').onclick = () => {
      inPreviewMode = true;
      currentAnswerPreview = 0;
      resultContainer.style.display = 'none';
      backToScoreBtn.style.display = 'inline-block';
      showAnswerPreview(currentAnswerPreview);
    };
  }
}


/*************************************************
 * 🔟 ANSWER PREVIEW MODE
 *************************************************/
function showAnswerPreview(index) {
  const q = quizData[index];

  quizContainer.innerHTML = `
    <h3>Question ${index + 1}</h3>
    <p>${q.question}</p>

    ${q.options.map((opt, i) => {
      let cls = '';
      if (i === q.answer) cls = 'correct';
      if (i === q.userAnswer && q.userAnswer !== q.answer) cls = 'incorrect';
      return `<div class="${cls}">${opt}</div>`;
    }).join('')}

    <div class="explanation">${q.explanation || ''}</div>
  `;

  prevBtn.style.display = index > 0 ? 'inline-block' : 'none';
  nextBtn.style.display = index < quizData.length - 1 ? 'inline-block' : 'none';
}


/*************************************************
 * 1️⃣1️⃣ BACK TO SCORE
 *************************************************/
backToScoreBtn.onclick = () => {
  quizContainer.innerHTML = '';
  resultContainer.style.display = 'block';
  backToScoreBtn.style.display = 'none';
  prevBtn.style.display = 'none';
  nextBtn.style.display = 'none';
};


/*************************************************
 * 1️⃣2️⃣ START EXAM INIT
 *************************************************/
startExam();

