/*************************************************
 * 1️⃣ SHUFFLE FUNCTION
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
 * 2️⃣ PICK RANDOM QUESTIONS & SHUFFLE OPTIONS
 *************************************************/
const TOTAL_QUESTIONS = 60;

// Function to generate fresh quiz data for new attempt
function generateQuizData() {
  let selected = shuffleArray(allQuizData).slice(0, TOTAL_QUESTIONS);
  return selected.map(q => {
    const correctText = q.options[q.answer];
    const newOptions = shuffleArray(q.options);
    const newAnswerIndex = newOptions.indexOf(correctText);
    return { ...q, options: newOptions, answer: newAnswerIndex, userAnswer: null };
  });
}

let quizData = generateQuizData();

/*************************************************
 * 3️⃣ CBT EXAM ENGINE ELEMENTS
 *************************************************/
const quizContainer = document.getElementById('quiz');
const resultContainer = document.getElementById('result');
const submitBtn = document.getElementById('submit');
const backToScoreBtn = document.getElementById('backToScore');
const prevBtn = document.getElementById('prev');
const nextBtn = document.getElementById('next');
let timerElement = document.getElementById('timer');

let timeLeft = 30 * 60;
let currentQuestion = 0;
let currentAnswerPreview = 0;
let inPreviewMode = false;
let timerInterval;
let quizStarted = false;

/*************************************************
 * 4️⃣ START SCREEN
 *************************************************/
function showStartScreen() {
  quizContainer.innerHTML = `
    <div id="introText" style="font-size:1.2rem;text-align:center;opacity:0;transition:opacity 1s ease;padding:40px;">
      <p id="typewriter"></p>
    </div>
  `;
  prevBtn.style.display = 'none';
  nextBtn.style.display = 'none';
  submitBtn.style.display = 'none';
  backToScoreBtn.style.display = 'none';
  resultContainer.innerHTML = '';
  timerElement.style.display = 'none';

  const text = "Get ready for Mock CBT Exam...";
  const el = document.getElementById('typewriter');
  let i = 0;
  const typing = setInterval(() => {
    el.textContent += text[i];
    i++;
    if (i === text.length) {
      clearInterval(typing);
      showStartButton();
    }
  }, 80);
  setTimeout(() => { document.getElementById('introText').style.opacity = '1'; }, 200);
}

function showStartButton() {
  const btn = document.createElement('button');
  btn.id = 'startExam';
  btn.textContent = 'Start Exam';
  btn.style.cssText = "padding:12px 25px;margin-top:20px;font-weight:bold;background:#1E90FF;color:#fff;border:none;border-radius:6px;cursor:pointer;opacity:0;transition:opacity 1s ease;";
  quizContainer.appendChild(btn);
  setTimeout(() => { btn.style.opacity = '1'; }, 500);
  btn.addEventListener('click', startExam);
}

/*************************************************
 * 5️⃣ START EXAM
 *************************************************/
function startExam() {
  quizStarted = true;
  currentQuestion = 0;
  timeLeft = 30 * 60;
  resultContainer.innerHTML = '';
  timerElement.style.display = 'block';
  showQuestion(currentQuestion);
  startTimer();
}

/*************************************************
 * 6️⃣ SHOW QUESTION
 *************************************************/
function showQuestion(index) {
  const q = quizData[index];
  const answers = q.options.map((opt, i) =>
    `<label style="display:block;text-align:left;margin:6px 0;">
      <input type='radio' name='question${index}' value='${i}' ${q.userAnswer == i ? 'checked' : ''}> ${opt}
    </label>`
  ).join('');
  quizContainer.innerHTML = `
    <div class='question'>
      <p><b>${index + 1}. ${q.question}</b></p>
      <div class='answers'>${answers}</div>
    </div>`;
  prevBtn.style.display = index === 0 ? 'none' : 'inline-block';
  nextBtn.style.display = index < quizData.length - 1 ? 'inline-block' : 'none';
  submitBtn.style.display = index === quizData.length - 1 ? 'inline-block' : 'none';
}

/*************************************************
 * 7️⃣ SAVE ANSWER
 *************************************************/
function saveAnswer() {
  const selected = document.querySelector(`input[name=question${currentQuestion}]:checked`);
  quizData[currentQuestion].userAnswer = selected ? parseInt(selected.value) : null;
}

function nextQuestion() {
  saveAnswer();
  if (currentQuestion < quizData.length - 1) {
    currentQuestion++;
    showQuestion(currentQuestion);
  }
}

function prevQuestion() {
  saveAnswer();
  if (currentQuestion > 0) {
    currentQuestion--;
    showQuestion(currentQuestion);
  }
}

/*************************************************
 * 8️⃣ SHOW RESULTS
 *************************************************/
function showResults() {
  saveAnswer();
  clearInterval(timerInterval);
  timerElement.style.display = 'none';

  let score = 0;
  quizData.forEach(q => { if (q.userAnswer == q.answer) score++; });

  const percent = (score / quizData.length) * 100;
  const grade =
    percent >= 90 ? 'A' :
    percent >= 70 ? 'B' :
    percent >= 50 ? 'C' : 'F';

  const gradeColor = grade === 'F' ? 'red' : 'inherit';
// 🔴 Show WhatsApp live tutorial ONLY for F grade
if (grade === 'A', 'B', 'C', 'F') {
  document.getElementById('waHelpItem').style.display = 'block';
} else {
  document.getElementById('waHelpItem').style.display = 'none';
}
  resultContainer.innerHTML = `
    <div style="font-size:1.1rem;">
      You scored <b>${score}</b> out of <b>${quizData.length}</b> (${percent.toFixed(0)}%)<br>
      Grade: <b style="color:${gradeColor}">${grade}</b><br><br>

      ${grade !== 'F'
        ? `<button id="showAnswersBtn">Preview Correct Answers & Explanation</button>`
        : ``}<br><br>

      ${grade === 'F'
        ? `<button id="retryBtn">Retry Exam</button>`
        : ``}

      <button id="homeBtn">Home</button>
    </div>
  `;

  quizContainer.innerHTML = '';
  prevBtn.style.display = 'none';
  nextBtn.style.display = 'none';
  submitBtn.style.display = 'none';
  backToScoreBtn.style.display = 'none';

  // Preview answers (only if not F)
  if (grade !== 'F') {
    document.getElementById('showAnswersBtn').addEventListener('click', () => {
      inPreviewMode = true;
      currentAnswerPreview = 0;
      showAnswerPreview(currentAnswerPreview);

      backToScoreBtn.style.display = 'inline-block';
      resultContainer.style.display = 'none';
    });
  }

  // Retry Exam (only if F)
  if (grade === 'F') {
    document.getElementById('retryBtn').addEventListener('click', () => {
      quizData = generateQuizData();
      showStartScreen();
    });
  }

  // Home button
  document.getElementById('homeBtn').addEventListener('click', () => {
    window.location.href = "index.html";
  });
}

/*************************************************
 * 9️⃣ SHOW ANSWER + EXPLANATION
 *************************************************/
function showAnswerPreview(index) {
  const q = quizData[index];
  const userAns = q.userAnswer != null ? q.options[q.userAnswer] : 'No Answer';
  quizContainer.innerHTML = `
    <div class='question'>
      <p><b>${index + 1}. ${q.question}</b></p>
      <p>Your Answer: <span class='${q.userAnswer == q.answer ? 'correct' : 'incorrect'}'>${userAns}</span></p>
      <p>Correct Answer: <span class='correct'>${q.options[q.answer]}</span></p>
      <p class='explanation'><b>Explanation:</b> ${q.explanation || 'No explanation provided.'}</p>
    </div>
  `;
  prevBtn.style.display = index === 0 ? 'none' : 'inline-block';
  nextBtn.style.display = index === quizData.length - 1 ? 'none' : 'inline-block';
  submitBtn.style.display = 'none';
}

/*************************************************
 * 10️⃣ NAVIGATION
 *************************************************/
prevBtn.addEventListener('click', () => {
  if (inPreviewMode) {
    if (currentAnswerPreview > 0) currentAnswerPreview--;
    showAnswerPreview(currentAnswerPreview);
  } else prevQuestion();
});

nextBtn.addEventListener('click', () => {
  if (inPreviewMode) {
    if (currentAnswerPreview < quizData.length - 1) currentAnswerPreview++;
    showAnswerPreview(currentAnswerPreview);
  } else nextQuestion();
});

submitBtn.addEventListener('click', showResults);

backToScoreBtn.addEventListener('click', () => {
  inPreviewMode = false;
  backToScoreBtn.style.display = 'none';
  quizContainer.innerHTML = '';
  resultContainer.style.display = 'block';
  prevBtn.style.display = 'none';
  nextBtn.style.display = 'none';
});

/*************************************************
 * 11️⃣ TIMER
 *************************************************/
function startTimer() {
  clearInterval(timerInterval);
  timerInterval = setInterval(() => {
    let minutes = Math.floor(timeLeft / 60);
    let seconds = timeLeft % 60;
    timerElement.textContent = `Time Left: ${minutes}:${seconds < 10 ? '0' : ''}${seconds}`;
    timeLeft--;
    if (timeLeft < 0) {
      clearInterval(timerInterval);
      showResults();
      alert('Time is up! Your answers have been submitted automatically.');
    }
  }, 1000);
}
const waModal = document.getElementById('waModal');
const whatsappLink = document.getElementById('whatsappLink');
const tryBtn = document.getElementById('tryBtn');
const closeBtn = document.getElementById('closeBtn');
const waNumber = document.getElementById('waNumber').textContent;

// Open modal
whatsappLink?.addEventListener('click', (e) => {
  e.preventDefault();
  waModal.style.display = 'flex';

  // Copy number automatically (offline-safe)
  navigator.clipboard?.writeText(waNumber).catch(() => {});
});

// Try opening WhatsApp (API – may fail offline)
tryBtn?.addEventListener('click', () => {
  const url = "https://wa.me/message/UKRVOLM35THWA1";
  window.open(url, "_self");
});

// Close modal
closeBtn?.addEventListener('click', () => {
  waModal.style.display = 'none';
});
