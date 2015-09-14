var gulp = require('gulp');
var rename = require('gulp-rename');
var uglify = require('gulp-uglify');

gulp.task('duplicate', function () {
  return gulp.src('./node_modules/dalliance/build/*.js')
    .pipe(gulp.dest('./dist'));
});

gulp.task('duplicate-css', function(){
  return gulp.src('./node_modules/dalliance/css/*')
    .pipe(gulp.dest('./css'));
});

gulp.task('duplicate-fonts', function(){
  return gulp.src('./node_modules/dalliance/fonts/*')
    .pipe(gulp.dest('./fonts'));
});

gulp.task('duplicate-img', function(){
  return gulp.src('./node_modules/dalliance/img/*')
    .pipe(gulp.dest('./img'));
});

gulp.task('compress', ['duplicate'], function () {
  return gulp.src(['./dist/*.js', '!./dist/*.min.js'])
    .pipe(uglify())
    .pipe(rename({
      extname: '.min.js'
    }))
    .pipe(gulp.dest('dist'));
});

gulp.task('default', ['duplicate', 'duplicate-css','duplicate-fonts','duplicate-img','compress']);