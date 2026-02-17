/*************************************************
 * 1️⃣ ADD YOUR QUESTIONS HERE (MANUALLY)
 * ➤ Each question must have:
 *    - question: string
 *    - options: array of strings
 *    - answer: index of correct option (0-based)
 *    - explanation: string
 *    - subject: string (e.g., "Math", "Physics", "English")
 *************************************************/

let quizData1 = [
  {
    question: "What is 2 + 2?",
    options: ["3", "4", "5", "6"],
    answer: 1,
    explanation: "2 + 2 equals 4.",
    subject: "Math"
  },
  {
    question: "What is the square root of 16?",
    options: ["2", "4", "8", "16"],
    answer: 1,
    explanation: "The square root of 16 is 4.",
    subject: "Math"
  }
];

let quizData2 = [
  {
    question: "Water boils at what temperature?",
    options: ["50°C", "100°C", "200°C", "0°C"],
    answer: 1,
    explanation: "Water boils at 100°C under normal pressure.",
    subject: "Physics"
  }
];

let quizData3 = [
  {
    question: "Who wrote 'Romeo and Juliet'?",
    options: ["Shakespeare", "Dickens", "Hemingway", "Orwell"],
    answer: 0,
    explanation: "William Shakespeare wrote 'Romeo and Juliet'.",
    subject: "English"
  }
];

let quizData4 = [
  {
    question: "What gas do plants release during photosynthesis?",
    options: ["Oxygen", "Carbon Dioxide", "Nitrogen", "Hydrogen"],
    answer: 0,
    explanation: "Plants release Oxygen during photosynthesis.",
    subject: "Biology"
  }
];

let quizData5 = [
  {
    question: "What is the capital of France?",
    options: ["Paris", "London", "Berlin", "Rome"],
    answer: 0,
    explanation: "Paris is the capital of France.",
    subject: "Geography"
  }
];

/*************************************************
 * 2️⃣ AUTO-COLLECT ALL QUIZ DATA
 *************************************************/

const allQuizData = [
  ...quizData1,
  ...quizData2,
  ...quizData3,
  ...quizData4,
  ...quizData5
];