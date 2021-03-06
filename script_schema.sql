USE [master]
GO
/****** Object:  Database [zakrzews_a]    Script Date: 01.02.2017 11:45:04 ******/
CREATE DATABASE [zakrzews_a]
 CONTAINMENT = NONE
 ON  PRIMARY 
( NAME = N'zakrzews_a', FILENAME = N'C:\Program Files\Microsoft SQL Server\MSSQL12.MSSQLSERVER\MSSQL\DATA\zakrzews_a.mdf' , SIZE = 10240KB , MAXSIZE = 30720KB , FILEGROWTH = 2048KB )
 LOG ON 
( NAME = N'zakrzews_a_log', FILENAME = N'C:\Program Files\Microsoft SQL Server\MSSQL12.MSSQLSERVER\MSSQL\DATA\zakrzews_a.ldf' , SIZE = 10240KB , MAXSIZE = 30720KB , FILEGROWTH = 2048KB )
GO
ALTER DATABASE [zakrzews_a] SET COMPATIBILITY_LEVEL = 120
GO
IF (1 = FULLTEXTSERVICEPROPERTY('IsFullTextInstalled'))
begin
EXEC [zakrzews_a].[dbo].[sp_fulltext_database] @action = 'enable'
end
GO
ALTER DATABASE [zakrzews_a] SET ANSI_NULL_DEFAULT OFF 
GO
ALTER DATABASE [zakrzews_a] SET ANSI_NULLS OFF 
GO
ALTER DATABASE [zakrzews_a] SET ANSI_PADDING OFF 
GO
ALTER DATABASE [zakrzews_a] SET ANSI_WARNINGS OFF 
GO
ALTER DATABASE [zakrzews_a] SET ARITHABORT OFF 
GO
ALTER DATABASE [zakrzews_a] SET AUTO_CLOSE OFF 
GO
ALTER DATABASE [zakrzews_a] SET AUTO_SHRINK OFF 
GO
ALTER DATABASE [zakrzews_a] SET AUTO_UPDATE_STATISTICS ON 
GO
ALTER DATABASE [zakrzews_a] SET CURSOR_CLOSE_ON_COMMIT OFF 
GO
ALTER DATABASE [zakrzews_a] SET CURSOR_DEFAULT  GLOBAL 
GO
ALTER DATABASE [zakrzews_a] SET CONCAT_NULL_YIELDS_NULL OFF 
GO
ALTER DATABASE [zakrzews_a] SET NUMERIC_ROUNDABORT OFF 
GO
ALTER DATABASE [zakrzews_a] SET QUOTED_IDENTIFIER OFF 
GO
ALTER DATABASE [zakrzews_a] SET RECURSIVE_TRIGGERS OFF 
GO
ALTER DATABASE [zakrzews_a] SET  ENABLE_BROKER 
GO
ALTER DATABASE [zakrzews_a] SET AUTO_UPDATE_STATISTICS_ASYNC OFF 
GO
ALTER DATABASE [zakrzews_a] SET DATE_CORRELATION_OPTIMIZATION OFF 
GO
ALTER DATABASE [zakrzews_a] SET TRUSTWORTHY OFF 
GO
ALTER DATABASE [zakrzews_a] SET ALLOW_SNAPSHOT_ISOLATION OFF 
GO
ALTER DATABASE [zakrzews_a] SET PARAMETERIZATION SIMPLE 
GO
ALTER DATABASE [zakrzews_a] SET READ_COMMITTED_SNAPSHOT OFF 
GO
ALTER DATABASE [zakrzews_a] SET HONOR_BROKER_PRIORITY OFF 
GO
ALTER DATABASE [zakrzews_a] SET RECOVERY SIMPLE 
GO
ALTER DATABASE [zakrzews_a] SET  MULTI_USER 
GO
ALTER DATABASE [zakrzews_a] SET PAGE_VERIFY CHECKSUM  
GO
ALTER DATABASE [zakrzews_a] SET DB_CHAINING OFF 
GO
ALTER DATABASE [zakrzews_a] SET FILESTREAM( NON_TRANSACTED_ACCESS = OFF ) 
GO
ALTER DATABASE [zakrzews_a] SET TARGET_RECOVERY_TIME = 0 SECONDS 
GO
ALTER DATABASE [zakrzews_a] SET DELAYED_DURABILITY = DISABLED 
GO
EXEC sys.sp_db_vardecimal_storage_format N'zakrzews_a', N'ON'
GO
USE [zakrzews_a]
GO
/****** Object:  UserDefinedFunction [dbo].[freePlacesForDay]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE FUNCTION [dbo].[freePlacesForDay] (@DayID int)
RETURNS int
AS
BEGIN
		DECLARE @allPlaces int;
		SET @allPlaces = (SELECT NumberOfPlaces FROM Conference WHERE ConfID = (SELECT ConfID FROM Days WHERE DayID = @DayID))

		DECLARE @takenPlaces int;
		SET @takenPlaces = (SELECT SUM(NumberOfPlaces) FROM Day_Registration WHERE DayID = @DayID)

		RETURN (@allPlaces - @takenPlaces)
END


GO
/****** Object:  UserDefinedFunction [dbo].[freePlacesForDayRegistration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
-- =============================================
-- Author:		<Author,,Name>
-- Create date: <Create Date, ,>
-- Description:	<Description, ,>
-- =============================================
CREATE FUNCTION [dbo].[freePlacesForDayRegistration] (@DayRegID int)
RETURNS int
AS
BEGIN
		DECLARE @allPlaces int;
		SET @allPlaces = (SELECT NumberOfPlaces FROM Day_Registration WHERE DayRegID = @DayRegID)

		DECLARE @takenPlaces int;
		SET @takenPlaces = (SELECT COUNT(*) FROM Which_Day WHERE DayRegID = @DayRegID)

		RETURN (@allPlaces - @takenPlaces)
END


GO
/****** Object:  UserDefinedFunction [dbo].[freePlacesForWorkshop]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
------------------------------
------------FUNKCJE-----------
------------------------------

--freePlacesForWorkshop
CREATE FUNCTION [dbo].[freePlacesForWorkshop] (@LessonID int)
RETURNS int
AS
BEGIN
		DECLARE @allPlaces int;
		SET @allPlaces = (SELECT NumberOfPlaces FROM Workshop_Lesson WHERE LessonID = @LessonID)

		DECLARE @takenPlaces int;
		SET @takenPlaces = (SELECT SUM(NumberOfPlaces) FROM Workshop_Registration WHERE LessonID = @LessonID)

		RETURN (@allPlaces - @takenPlaces)
END


GO
/****** Object:  UserDefinedFunction [dbo].[freePlacesForWorkshopRegistration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO


--freePlacesForWorkshopRegistration
CREATE FUNCTION [dbo].[freePlacesForWorkshopRegistration] (@WorkRegID int)
RETURNS int
AS
BEGIN
		DECLARE @allPlaces int;
		SET @allPlaces = (SELECT NumberOfPlaces FROM Workshop_Registration WHERE WorkRegID = @WorkRegID)

		DECLARE @takenPlaces int;
		SET @takenPlaces = (SELECT COUNT(*) FROM Which_Workshop WHERE WorkRegID = @WorkRegID)

		RETURN (@allPlaces - @takenPlaces)
END


GO
/****** Object:  UserDefinedFunction [dbo].[numberOfPlacesNotWorkshoped]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE FUNCTION [dbo].[numberOfPlacesNotWorkshoped] (@DayID int)
RETURNS int
AS
BEGIN
		DECLARE @allPlaces int;
		SET @allPlaces = (	SELECT NumberOfPlaces 
							FROM Conference 
							WHERE ConfID = (SELECT ConfID 
											FROM Days
											WHERE DayID = @DayID));

		DECLARE @takenPlaces int;
		SET @takenPlaces = (SELECT SUM(NumberOfPlaces)
							FROM Workshop_Lesson WL
							WHERE DayID = @DayID);

		RETURN (@allPlaces - ISNULL(@takenPlaces, 0))
END

GO
/****** Object:  UserDefinedFunction [dbo].[numberOfStudentsInRegistration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE FUNCTION [dbo].[numberOfStudentsInRegistration](@RegistrationID int)
RETURNS int
AS
BEGIN

		DECLARE @StudentsNumber int = (	SELECT COUNT(*)
										FROM Student S
										INNER JOIN Participant P ON P.ParticipantID = S.ParticipantID
										INNER JOIN Which_Day WD ON WD.ParticipantID = P.ParticipantID
										INNER JOIN Day_Registration DR ON DR.DayRegID = WD.DayRegID
										WHERE RegistrationID = @RegistrationID)

		RETURN (@StudentsNumber)
END
GO
/****** Object:  UserDefinedFunction [dbo].[registrationCost]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO

CREATE FUNCTION [dbo].[registrationCost] (@RegistrationID int)
RETURNS float
AS
BEGIN

	DECLARE @regDate date = (
		SELECT RegistrationDate from Registration
		WHERE RegistrationID = @RegistrationID
	)

	DECLARE @discount float = (
		SELECT D.Discount FROM Discount as D
		WHERE D.StartDate <= @regDate AND D.EndDate >= @regDate AND D.ConfID IN (
			SELECT Da.ConfID FROM Days as Da
			WHERE Da.DayID IN (
				SELECT DR.DayID FROM Day_Registration as DR
				WHERE DR.RegistrationID = @RegistrationID
			)
		)
	)
	IF(@discount IS NULL)
		BEGIN
		;SET @discount = 0
		END

	DECLARE @dayCost float =(1-@discount)*(
		SELECT SUM(PRICE) FROM(
			SELECT /*DR.DayID as ID, */((DR.NumberOfPlaces-dbo.numberOfStudentsInRegistration(@RegistrationID))*D.DayPrice
			+dbo.numberOfStudentsInRegistration(@RegistrationID)*(1-C.StudentDiscount)*D.DayPrice) as PRICE FROM Day_Registration as DR
			INNER JOIN Days as D
			ON D.DayID = DR.DayID
			INNER JOIN Conference as C
			ON C.ConfID = D.ConfID
			WHERE DR.RegistrationID = @RegistrationID
		) AS ID
		--GROUP BY PRICE
	)

	DECLARE @workshopCost float = (
		SELECT SUM(PRICE) FROM(
			SELECT WR.NumberOfPlaces * WL.Price as PRICE FROM Workshop_Registration as WR
			INNER JOIN Workshop_Lesson as WL
			ON WL.LessonID = WR.LessonID
			INNER JOIN Day_Registration as DR
			ON DR.DayRegID = WR.DayRegID
			WHERE DR.RegistrationID = @RegistrationID
		) AS ID
		--GROUP BY PRICE
	)

	IF(@workshopCost IS NULL)
		BEGIN
		;SET @workshopCost = 0
		END
	IF(@dayCost IS NULL)
		BEGIN
		;SET @dayCost = 0
		END
				
	RETURN(@dayCost + @workshopCost)
	
END
GO
/****** Object:  UserDefinedFunction [dbo].[workshopCollision]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO

--workshop collision
CREATE FUNCTION [dbo].[workshopCollision] (@lessonID1 int, @lessonID2 int)
RETURNS bit
AS
BEGIN
		DECLARE @startTime1 datetime;
		DECLARE @endTime1 datetime;
		DECLARE @startTime2 datetime;
		DECLARE @endTime2 datetime;

		SET @startTime1 = (SELECT startHour FROM Workshop_Lesson WHERE LessonID =@lessonID1)
		SET @endTime1 = (SELECT EndHour FROM Workshop_Lesson WHERE LessonID =@lessonID1)
		SET @startTime2 = (SELECT startHour FROM Workshop_Lesson WHERE LessonID =@lessonID2)
		SET @endTime2 = (SELECT EndHour FROM Workshop_Lesson WHERE LessonID =@lessonID2)

		DECLARE @collision bit;
		IF((@startTime1 > @startTime2 AND @startTime1 < @endTime2) OR (@endTime1 > @startTime2 AND @endTime1 < @endTime2))
			SET @collision = 1
		ELSE
			SET @collision = 0
		
		RETURN @collision 

END

GO
/****** Object:  Table [dbo].[Client]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
SET ANSI_PADDING ON
GO
CREATE TABLE [dbo].[Client](
	[ClientID] [int] IDENTITY(1,1) NOT NULL,
	[FirstName] [varchar](50) NOT NULL,
	[LastName] [varchar](50) NOT NULL,
	[Phone] [varchar](15) NOT NULL,
	[E-mail] [varchar](255) NOT NULL,
 CONSTRAINT [PK_Client] PRIMARY KEY CLUSTERED 
(
	[ClientID] ASC
)WITH (PAD_INDEX = OFF, STATISTICS_NORECOMPUTE = OFF, IGNORE_DUP_KEY = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY]
) ON [PRIMARY]

GO
SET ANSI_PADDING OFF
GO
/****** Object:  Table [dbo].[Company]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
SET ANSI_PADDING ON
GO
CREATE TABLE [dbo].[Company](
	[CompanyID] [int] IDENTITY(1,1) NOT NULL,
	[ClientID] [int] NOT NULL,
	[CompanyName] [varchar](50) NOT NULL,
	[City] [varchar](50) NOT NULL,
	[Address] [varchar](50) NOT NULL,
	[Country] [varchar](50) NOT NULL,
 CONSTRAINT [PK_Company] PRIMARY KEY CLUSTERED 
(
	[CompanyID] ASC
)WITH (PAD_INDEX = OFF, STATISTICS_NORECOMPUTE = OFF, IGNORE_DUP_KEY = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY]
) ON [PRIMARY]

GO
SET ANSI_PADDING OFF
GO
/****** Object:  Table [dbo].[Conference]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
SET ANSI_PADDING ON
GO
CREATE TABLE [dbo].[Conference](
	[ConfID] [int] IDENTITY(1,1) NOT NULL,
	[Name] [varchar](255) NOT NULL,
	[StartDate] [date] NOT NULL,
	[EndDate] [date] NOT NULL,
	[City] [varchar](50) NOT NULL,
	[Address] [varchar](50) NOT NULL,
	[Country] [varchar](50) NOT NULL,
	[NumberOfPlaces] [int] NOT NULL,
	[StudentDiscount] [float] NOT NULL,
 CONSTRAINT [PK_Conference] PRIMARY KEY CLUSTERED 
(
	[ConfID] ASC
)WITH (PAD_INDEX = OFF, STATISTICS_NORECOMPUTE = OFF, IGNORE_DUP_KEY = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY]
) ON [PRIMARY]

GO
SET ANSI_PADDING OFF
GO
/****** Object:  Table [dbo].[Day_Registration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE TABLE [dbo].[Day_Registration](
	[DayRegID] [int] IDENTITY(1,1) NOT NULL,
	[RegistrationID] [int] NOT NULL,
	[DayID] [int] NOT NULL,
	[NumberOfPlaces] [int] NOT NULL,
 CONSTRAINT [PK_Day_Registration] PRIMARY KEY CLUSTERED 
(
	[DayRegID] ASC
)WITH (PAD_INDEX = OFF, STATISTICS_NORECOMPUTE = OFF, IGNORE_DUP_KEY = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY]
) ON [PRIMARY]

GO
/****** Object:  Table [dbo].[Days]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE TABLE [dbo].[Days](
	[DayID] [int] IDENTITY(1,1) NOT NULL,
	[ConfID] [int] NOT NULL,
	[Date] [date] NOT NULL,
	[DayPrice] [money] NOT NULL,
 CONSTRAINT [PK_Days] PRIMARY KEY CLUSTERED 
(
	[DayID] ASC
)WITH (PAD_INDEX = OFF, STATISTICS_NORECOMPUTE = OFF, IGNORE_DUP_KEY = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY]
) ON [PRIMARY]

GO
/****** Object:  Table [dbo].[Discount]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE TABLE [dbo].[Discount](
	[DiscountID] [int] IDENTITY(1,1) NOT NULL,
	[ConfID] [int] NOT NULL,
	[Discount] [float] NOT NULL,
	[StartDate] [date] NOT NULL,
	[EndDate] [date] NOT NULL,
 CONSTRAINT [PK_Discount] PRIMARY KEY CLUSTERED 
(
	[DiscountID] ASC
)WITH (PAD_INDEX = OFF, STATISTICS_NORECOMPUTE = OFF, IGNORE_DUP_KEY = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY]
) ON [PRIMARY]

GO
/****** Object:  Table [dbo].[Participant]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
SET ANSI_PADDING ON
GO
CREATE TABLE [dbo].[Participant](
	[ParticipantID] [int] IDENTITY(1,1) NOT NULL,
	[FirstName] [varchar](50) NULL,
	[LastName] [varchar](50) NULL,
 CONSTRAINT [PK_Participant] PRIMARY KEY CLUSTERED 
(
	[ParticipantID] ASC
)WITH (PAD_INDEX = OFF, STATISTICS_NORECOMPUTE = OFF, IGNORE_DUP_KEY = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY]
) ON [PRIMARY]

GO
SET ANSI_PADDING OFF
GO
/****** Object:  Table [dbo].[Registration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE TABLE [dbo].[Registration](
	[RegistrationID] [int] IDENTITY(1,1) NOT NULL,
	[ConfID] [int] NOT NULL,
	[ClientID] [int] NOT NULL,
	[RegistrationDate] [date] NOT NULL,
	[PayDate] [date] NULL,
 CONSTRAINT [PK_Registration] PRIMARY KEY CLUSTERED 
(
	[RegistrationID] ASC
)WITH (PAD_INDEX = OFF, STATISTICS_NORECOMPUTE = OFF, IGNORE_DUP_KEY = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY]
) ON [PRIMARY]

GO
/****** Object:  Table [dbo].[Student]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
SET ANSI_PADDING ON
GO
CREATE TABLE [dbo].[Student](
	[ParticipantID] [int] NOT NULL,
	[StudentCardNumber] [varchar](50) NOT NULL
) ON [PRIMARY]

GO
SET ANSI_PADDING OFF
GO
/****** Object:  Table [dbo].[Which_Day]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE TABLE [dbo].[Which_Day](
	[ParticipantID] [int] NOT NULL,
	[DayRegID] [int] NOT NULL
) ON [PRIMARY]

GO
/****** Object:  Table [dbo].[Which_Workshop]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE TABLE [dbo].[Which_Workshop](
	[ParticipantID] [int] NOT NULL,
	[WorkRegID] [int] NOT NULL
) ON [PRIMARY]

GO
/****** Object:  Table [dbo].[Workshop_Lesson]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE TABLE [dbo].[Workshop_Lesson](
	[LessonID] [int] IDENTITY(1,1) NOT NULL,
	[WorkshopID] [int] NOT NULL,
	[DayID] [int] NOT NULL,
	[NumberOfPlaces] [int] NOT NULL,
	[Price] [money] NOT NULL,
	[StartHour] [datetime] NOT NULL,
	[EndHour] [datetime] NOT NULL,
 CONSTRAINT [PK_Workshop_Lesson] PRIMARY KEY CLUSTERED 
(
	[LessonID] ASC
)WITH (PAD_INDEX = OFF, STATISTICS_NORECOMPUTE = OFF, IGNORE_DUP_KEY = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY]
) ON [PRIMARY]

GO
/****** Object:  Table [dbo].[Workshop_Name]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
SET ANSI_PADDING ON
GO
CREATE TABLE [dbo].[Workshop_Name](
	[WorkshopID] [int] IDENTITY(1,1) NOT NULL,
	[Title] [varchar](50) NOT NULL,
	[Description] [varchar](255) NOT NULL,
 CONSTRAINT [PK_Workshop_Name] PRIMARY KEY CLUSTERED 
(
	[WorkshopID] ASC
)WITH (PAD_INDEX = OFF, STATISTICS_NORECOMPUTE = OFF, IGNORE_DUP_KEY = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY]
) ON [PRIMARY]

GO
SET ANSI_PADDING OFF
GO
/****** Object:  Table [dbo].[Workshop_Registration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE TABLE [dbo].[Workshop_Registration](
	[WorkRegID] [int] IDENTITY(1,1) NOT NULL,
	[DayRegID] [int] NOT NULL,
	[LessonID] [int] NOT NULL,
	[Date] [date] NOT NULL,
	[NumberOfPlaces] [int] NOT NULL,
 CONSTRAINT [PK_Workshop_Registration] PRIMARY KEY CLUSTERED 
(
	[WorkRegID] ASC
)WITH (PAD_INDEX = OFF, STATISTICS_NORECOMPUTE = OFF, IGNORE_DUP_KEY = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY]
) ON [PRIMARY]

GO
/****** Object:  View [dbo].[VIEW_ClientsRegistrations]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
/*ILOSC REJESTRACJI DANEGO KLIENTA--*/
CREATE VIEW [dbo].[VIEW_ClientsRegistrations]
AS
SELECT        C.ClientID, COUNT(*) AS NUMBER
FROM            dbo.Client AS C INNER JOIN
                         dbo.Registration AS R ON R.ClientID = C.ClientID
GROUP BY C.ClientID

GO
/****** Object:  View [dbo].[VIEW_ClientsToCall]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE VIEW [dbo].[VIEW_ClientsToCall]
AS
SELECT        C.ClientID, C.FirstName, C.LastName, C.phone, 'DAY' AS registrationType, R.RegistrationID, dbo.freePlacesForDayRegistration(DR.DayRegID) AS freePlaces
FROM            Client C INNER JOIN
                         Registration R ON R.ClientID = C.ClientID INNER JOIN
                         Day_Registration DR ON DR.RegistrationID = R.RegistrationID INNER JOIN
                         Days D ON D .DayID = DR.DayID INNER JOIN
                         Conference CO ON CO.ConfID = D .ConfID
WHERE        dbo.freePlacesForDayRegistration(DR.DayID) > 0 AND DATEDIFF(day, CONVERT(date, GETDATE()), CO.StartDate) BETWEEN 0 AND 14 
UNION
SELECT        C.ClientID, C.FirstName, C.LastName, C.phone, 'WORKSHOP' AS registrationType, R.RegistrationID, dbo.freePlacesForWorkshopRegistration(WR.WorkRegID) AS freePlaces
FROM            Client C INNER JOIN
                         Registration R ON R.ClientID = C.ClientID INNER JOIN
                         Day_Registration DR ON DR.RegistrationID = R.RegistrationID INNER JOIN
                         Days D ON D .DayID = DR.DayID INNER JOIN
                         Conference CO ON CO.ConfID = D .ConfID INNER JOIN
                         Workshop_Registration WR ON WR.DayRegID = DR.DayRegID
WHERE        dbo.freePlacesForWorkshopRegistration(WR.WorkRegID) > 0 AND DATEDIFF(day, CONVERT(date, GETDATE()), CO.StartDate) BETWEEN 0 AND 14 

GO
/****** Object:  View [dbo].[VIEW_CustomersThatShouldPayTomorrow]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO

create view [dbo].[VIEW_CustomersThatShouldPayTomorrow] as
select C.ClientID, C.FirstName, C.LastName, C.[E-mail], C.Phone
	from Client as C
inner join Registration as R
	on C.ClientID = R.ClientID
where DATEDIFF(day, convert(date, getdate()), R.RegistrationDate) = 6

GO
/****** Object:  View [dbo].[VIEW_OutdatedPayments]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO

CREATE VIEW [dbo].[VIEW_OutdatedPayments]
AS
SELECT RegistrationID FROM Registration
WHERE DATEDIFF(day, CONVERT(date, GETDATE()), RegistrationDate) > 7
GO
/****** Object:  View [dbo].[VIEW_ParticipantStats]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO


--statystyki uczestnika na ilu byl warsztatach i na ilu byl dniach

CREATE VIEW [dbo].[VIEW_ParticipantStats] AS

SELECT P.ParticipantID, COUNT(WD.DayRegID) AS number, 'DAYS' AS activityType
FROM Participant P
INNER JOIN Which_Day WD ON WD.ParticipantID = P.ParticipantID
GROUP BY P.ParticipantID

UNION

SELECT P.ParticipantID, COUNT(WW.WorkRegID) AS number, 'WORKSHOP' AS activityType
FROM Participant P
INNER JOIN Which_Workshop WW ON WW.ParticipantID = P.ParticipantID
GROUP BY P.ParticipantID

GO
/****** Object:  View [dbo].[VIEW_studentParticipants]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO

CREATE VIEW [dbo].[VIEW_studentParticipants]
AS
SELECT ParticipantID FROM Student
GO
ALTER TABLE [dbo].[Company]  WITH CHECK ADD  CONSTRAINT [FK_Company_Client] FOREIGN KEY([ClientID])
REFERENCES [dbo].[Client] ([ClientID])
GO
ALTER TABLE [dbo].[Company] CHECK CONSTRAINT [FK_Company_Client]
GO
ALTER TABLE [dbo].[Day_Registration]  WITH CHECK ADD  CONSTRAINT [FK_Day_Registration_Days] FOREIGN KEY([DayID])
REFERENCES [dbo].[Days] ([DayID])
GO
ALTER TABLE [dbo].[Day_Registration] CHECK CONSTRAINT [FK_Day_Registration_Days]
GO
ALTER TABLE [dbo].[Day_Registration]  WITH CHECK ADD  CONSTRAINT [FK_Day_Registration_Registration] FOREIGN KEY([RegistrationID])
REFERENCES [dbo].[Registration] ([RegistrationID])
GO
ALTER TABLE [dbo].[Day_Registration] CHECK CONSTRAINT [FK_Day_Registration_Registration]
GO
ALTER TABLE [dbo].[Days]  WITH CHECK ADD  CONSTRAINT [FK_Days_Conference] FOREIGN KEY([ConfID])
REFERENCES [dbo].[Conference] ([ConfID])
GO
ALTER TABLE [dbo].[Days] CHECK CONSTRAINT [FK_Days_Conference]
GO
ALTER TABLE [dbo].[Discount]  WITH CHECK ADD  CONSTRAINT [FK_Discount_Conference1] FOREIGN KEY([ConfID])
REFERENCES [dbo].[Conference] ([ConfID])
GO
ALTER TABLE [dbo].[Discount] CHECK CONSTRAINT [FK_Discount_Conference1]
GO
ALTER TABLE [dbo].[Registration]  WITH CHECK ADD  CONSTRAINT [FK_Registration_Client] FOREIGN KEY([ClientID])
REFERENCES [dbo].[Client] ([ClientID])
GO
ALTER TABLE [dbo].[Registration] CHECK CONSTRAINT [FK_Registration_Client]
GO
ALTER TABLE [dbo].[Registration]  WITH CHECK ADD  CONSTRAINT [FK_Registration_Conference] FOREIGN KEY([ConfID])
REFERENCES [dbo].[Conference] ([ConfID])
GO
ALTER TABLE [dbo].[Registration] CHECK CONSTRAINT [FK_Registration_Conference]
GO
ALTER TABLE [dbo].[Student]  WITH CHECK ADD  CONSTRAINT [FK_Student_Participant] FOREIGN KEY([ParticipantID])
REFERENCES [dbo].[Participant] ([ParticipantID])
GO
ALTER TABLE [dbo].[Student] CHECK CONSTRAINT [FK_Student_Participant]
GO
ALTER TABLE [dbo].[Which_Day]  WITH CHECK ADD  CONSTRAINT [FK_Which_Day_Day_Registration] FOREIGN KEY([DayRegID])
REFERENCES [dbo].[Day_Registration] ([DayRegID])
GO
ALTER TABLE [dbo].[Which_Day] CHECK CONSTRAINT [FK_Which_Day_Day_Registration]
GO
ALTER TABLE [dbo].[Which_Day]  WITH CHECK ADD  CONSTRAINT [FK_Which_Day_Participant] FOREIGN KEY([ParticipantID])
REFERENCES [dbo].[Participant] ([ParticipantID])
GO
ALTER TABLE [dbo].[Which_Day] CHECK CONSTRAINT [FK_Which_Day_Participant]
GO
ALTER TABLE [dbo].[Which_Workshop]  WITH CHECK ADD  CONSTRAINT [FK_Which_Workshop_Participant] FOREIGN KEY([ParticipantID])
REFERENCES [dbo].[Participant] ([ParticipantID])
GO
ALTER TABLE [dbo].[Which_Workshop] CHECK CONSTRAINT [FK_Which_Workshop_Participant]
GO
ALTER TABLE [dbo].[Which_Workshop]  WITH CHECK ADD  CONSTRAINT [FK_Which_Workshop_Workshop_Registration] FOREIGN KEY([WorkRegID])
REFERENCES [dbo].[Workshop_Registration] ([WorkRegID])
GO
ALTER TABLE [dbo].[Which_Workshop] CHECK CONSTRAINT [FK_Which_Workshop_Workshop_Registration]
GO
ALTER TABLE [dbo].[Workshop_Lesson]  WITH CHECK ADD  CONSTRAINT [FK_Workshop_Lesson_Days] FOREIGN KEY([DayID])
REFERENCES [dbo].[Days] ([DayID])
GO
ALTER TABLE [dbo].[Workshop_Lesson] CHECK CONSTRAINT [FK_Workshop_Lesson_Days]
GO
ALTER TABLE [dbo].[Workshop_Lesson]  WITH CHECK ADD  CONSTRAINT [FK_Workshop_Lesson_Workshop_Name] FOREIGN KEY([WorkshopID])
REFERENCES [dbo].[Workshop_Name] ([WorkshopID])
GO
ALTER TABLE [dbo].[Workshop_Lesson] CHECK CONSTRAINT [FK_Workshop_Lesson_Workshop_Name]
GO
ALTER TABLE [dbo].[Workshop_Registration]  WITH CHECK ADD  CONSTRAINT [FK_Workshop_Registration_Day_Registration] FOREIGN KEY([DayRegID])
REFERENCES [dbo].[Day_Registration] ([DayRegID])
GO
ALTER TABLE [dbo].[Workshop_Registration] CHECK CONSTRAINT [FK_Workshop_Registration_Day_Registration]
GO
ALTER TABLE [dbo].[Workshop_Registration]  WITH CHECK ADD  CONSTRAINT [FK_Workshop_Registration_Workshop_Lesson] FOREIGN KEY([LessonID])
REFERENCES [dbo].[Workshop_Lesson] ([LessonID])
GO
ALTER TABLE [dbo].[Workshop_Registration] CHECK CONSTRAINT [FK_Workshop_Registration_Workshop_Lesson]
GO
ALTER TABLE [dbo].[Conference]  WITH CHECK ADD  CONSTRAINT [CheckNumberOfPlaces_C] CHECK  (([NumberOfPlaces]>(0)))
GO
ALTER TABLE [dbo].[Conference] CHECK CONSTRAINT [CheckNumberOfPlaces_C]
GO
ALTER TABLE [dbo].[Conference]  WITH CHECK ADD  CONSTRAINT [CheckStartEnd_C] CHECK  (([StartDate]<=[EndDate]))
GO
ALTER TABLE [dbo].[Conference] CHECK CONSTRAINT [CheckStartEnd_C]
GO
ALTER TABLE [dbo].[Day_Registration]  WITH CHECK ADD  CONSTRAINT [CheckNumberOfPlaces_DR] CHECK  (([NumberOfPlaces]>(0)))
GO
ALTER TABLE [dbo].[Day_Registration] CHECK CONSTRAINT [CheckNumberOfPlaces_DR]
GO
ALTER TABLE [dbo].[Days]  WITH CHECK ADD  CONSTRAINT [CheckPrice] CHECK  (([DayPrice]>=(0)))
GO
ALTER TABLE [dbo].[Days] CHECK CONSTRAINT [CheckPrice]
GO
ALTER TABLE [dbo].[Discount]  WITH CHECK ADD  CONSTRAINT [CheckDiscount] CHECK  (([Discount]>=(0) AND [Discount]<=(1)))
GO
ALTER TABLE [dbo].[Discount] CHECK CONSTRAINT [CheckDiscount]
GO
ALTER TABLE [dbo].[Discount]  WITH CHECK ADD  CONSTRAINT [CheckStartEnd_Disc] CHECK  (([StartDate]<=[EndDate]))
GO
ALTER TABLE [dbo].[Discount] CHECK CONSTRAINT [CheckStartEnd_Disc]
GO
ALTER TABLE [dbo].[Workshop_Lesson]  WITH CHECK ADD  CONSTRAINT [CheckLessonHours] CHECK  (([StartHour]<[EndHour]))
GO
ALTER TABLE [dbo].[Workshop_Lesson] CHECK CONSTRAINT [CheckLessonHours]
GO
ALTER TABLE [dbo].[Workshop_Lesson]  WITH CHECK ADD  CONSTRAINT [CheckNumberOfPlaces] CHECK  (([NumberOfPlaces]>(0)))
GO
ALTER TABLE [dbo].[Workshop_Lesson] CHECK CONSTRAINT [CheckNumberOfPlaces]
GO
ALTER TABLE [dbo].[Workshop_Lesson]  WITH CHECK ADD  CONSTRAINT [CheckWorkshopPrice] CHECK  (([Price]>=(0)))
GO
ALTER TABLE [dbo].[Workshop_Lesson] CHECK CONSTRAINT [CheckWorkshopPrice]
GO
ALTER TABLE [dbo].[Workshop_Registration]  WITH CHECK ADD  CONSTRAINT [CheckNumberOfPlaces_WR] CHECK  (([NumberOfPlaces]>(0)))
GO
ALTER TABLE [dbo].[Workshop_Registration] CHECK CONSTRAINT [CheckNumberOfPlaces_WR]
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addClient]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_addClient]
@firstName varchar(50), @lastName varchar(50), @Phone varchar(15), @mail varchar(255)
AS
BEGIN
SET NOCOUNT ON;
		
		INSERT INTO Client(FirstName, LastName, Phone, [E-mail])
		VALUES(@firstName, @lastName, @Phone, @mail)
END

GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addCompany]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_addCompany]
@ClientID int, @CompanyName varchar(50), @City varchar(50), @address varchar(50), @country varchar(50)
AS
BEGIN
SET NOCOUNT ON;

		DECLARE @ClientID2 int = (SELECT ClientID FROM CLient WHERE CLientID = @ClientID)

		IF(@ClientID2 IS NULL)
		BEGIN
		;THROW 52000, 'There is no such client.',1
		END

		INSERT INTO Company(ClientID, CompanyName, City, Address, Country)
		VALUES(@ClientID, @CompanyName, @City, @address, @country)

END

GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addConference]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_addConference]
@Name varchar(255), @StartDate date, @EndDate date, @City varchar(50), @Address varchar(50), @Country varchar(50), @NumberOfPlaces int, @StudentDiscount float
AS
BEGIN
SET NOCOUNT ON;

		IF(@NumberOfPlaces < 0)
		BEGIN
		;THROW 52000, 'The number of places must not be negative.', 1
		END
		IF(@StartDate > @EndDate)
		BEGIN
		;THROW 52000, 'EndDate should not be earlier than StartDate.',1
		END
		IF(@StudentDiscount < 0 OR @StudentDiscount > 1)
		BEGIN
		;THROW 52000, 'The discount must be between 0 and 1.',1
		END

		INSERT INTO Conference(Name, StartDate, EndDate, City, Address, Country, studentDiscount, NumberOfPlaces)
		VALUES(@Name, @StartDate, @EndDate, @City, @Address, @Country, @StudentDiscount, @NumberOfPlaces)
END

GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addDay]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_addDay]
@ConfID int, @Date date, @DayPrice money
AS
BEGIN
SET NOCOUNT ON;
		
		DECLARE @Conference int = (SELECT ConfID FROM Conference WHERE ConfID = @ConfID)
		IF(@conference IS NULL)
		BEGIN
		;THROW 52000, 'There is no such conference.', 1
		END

		IF(@DayPrice < 0)
		BEGIN
		;THROW 52000, 'The price of day must not be negative.', 1
		END

		DECLARE @startDate date = (SELECT StartDate FROM Conference WHERE ConfID = @ConfID)
		DECLARE @EndDate date = (SELECT EndDate FROM Conference WHERE ConfID = @ConfID)
		IF(@Date < @startDate OR  @Date > @EndDate)
		BEGIN
		;THROW 52000, 'Date must be between startdate and enddate of conference.', 1
		END

		DECLARE @Date2 date = (SELECT date FROM Days WHERE ConfID = @ConfID AND Date = @Date)
		IF(@Date2 IS NOT NULL)
		BEGIN
		;THROW 52000, 'This day of conference already exists.', 1
		END

		INSERT INTO Days(ConfID, date, DayPrice)
		VALUES(@ConfID, @Date, @DayPrice)

END

GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addDayRegistration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_addDayRegistration]
@RegistrationID int, @NumberOfPlaces int, @DayID int
AS
BEGIN
SET NOCOUNT ON;
		
		DECLARE @Day int = (SELECT DayID FROM Days WHERE DayID = @DayID);
		IF(@Day IS NULL)
		BEGIN
		;THROW 52000, 'There is no such day.',1
		END

		DECLARE @Registration int = (SELECT RegistrationID FROM Registration WHERE RegistrationID = @RegistrationID);
		IF(@Registration IS NULL)
		BEGIN
		;THROW 52000, 'There is no such Registration.', 1
		END

		DECLARE @ConfID int = (SELECT ConfID FROM Registration WHERE RegistrationID = @RegistrationID)
		DECLARE @ConfID2 int = (SELECT COnfID FROM DAys WHERE DayID = @DayID)

		IF(@ConfID != @ConfID2)
		BEGIN
		;THROW 52000, 'You picked a day from wrong conference', 1
		END

		DECLARE @freePlacesforDay int = dbo.freePlacesForDay(@DayID);

		IF(@NumberOfPlaces > @freePlacesforDay)
		BEGIN
		PRINT 'Only ' + CAST(@freePlacesforDay AS VARCHAR(10)) + ' available'
		;THROW 52000, 'Not enouch places available', 1
		END

		INSERT INTO Day_Registration(RegistrationID, DayID, NumberOfPlaces)
		VALUES(@RegistrationID, @Day, @NumberOfPlaces)
END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addDiscount]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_addDiscount]
@ConfID int , @Discount float, @Startdate date, @EndDate date
AS
BEGIN
		
		DECLARE @Conf int =(SELECT COnfID FROM Conference WHERE ConfID = @ConfID)
		IF(@Conf IS NULL)
		BEGIN
		;THROW 52000, 'There is no such conference.', 1
		END

		IF(@Discount NOT BETWEEN 0 AND 1)
		BEGIN
		;THROW 52000, 'Discount should be between 0 and 1.', 1
		END

		IF(@StartDate > @EndDate)
		BEGIN
		;THROW 52000, 'StartDate should be before EndDate.', 1
		END


		IF(@ConfID  IN (	SELECT ConfID
							FROM Discount
							WHERE (@Startdate > StartDate AND @StartDate < EndDate) OR (@EndDate > StartDate AND @EndDate < EndDate)))
		BEGIN
		;THROW 52000, 'Collision with another discount from this conference', 1
		END

		INSERT INTO Discount(ConfID, Discount, StartDate, EndDate)
		VALUES(@ConfID, @Discount, @Startdate, @EndDate)

END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addLesson]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_addLesson]
@WorkshopID int, @DayID int, @NumberOfPlaces int, @Price money, @StartHour datetime, @EndHour datetime
AS
BEGIN
SET NOCOUNT ON;

		DECLARE @WorkshopID2 int = (SELECT WorkshopID FROM Workshop_Name WHERE WorkshopID = @WorkshopID)
		IF(@WorkshopID2 IS NULL)
		BEGIN
		;THROW 52000, 'There is no such workshop.', 1
		END

		DECLARE @DayID2 int = (SELECT DayID FROM Days WHERE DayID = @DayID)
		IF(@DayID2 IS NULL)
		BEGIN
		;THROW 52000, 'There is no such day.', 1
		END

		IF(@NumberOfPlaces < 0)
		BEGIN
		;THROW 52000, 'The number of places must not be negative.', 1
		END
		
		DECLARE @conferencePlaces int = (	SELECT C.numberOfPlaces
											FROM Conference C
											INNER JOIN Days D On D.ConfID = C.ConfID
											WHERE D.DayID = @DayID)

		IF(@conferencePlaces < @NumberOfPlaces)
		BEGIN
		PRINT 'Only ' + CAST(@conferencePlaces AS VARCHAR(10)) + ' Available'
		;THROW 52000, 'Not enough places available only', 1
		END

		IF(@StartHour > @EndHour)
		BEGIN
		;THROW 52000, 'EndHour should not be earlier than StartHour.', 1 
		END

		DECLARE @daydate date;
		SET @daydate = (SELECT date FROM Days WHERE DayID = @DayID)

		DECLARE @startDate date = CONVERT(date, @starthour)
		DECLARE @endDate date = CONVERT(date, @endhour)

		IF(@startDate != @daydate OR @endDate != @daydate)
		BEGIN
		;THROW 52000, 'You picked the wrong start hour or endhour - not from this day', 1
		END

		INSERT INTO Workshop_Lesson(WorkshopID, DayID, NumberOfPlaces, Price, StartHour, EndHour)
		VALUES(@WorkshopID, @DayID, @NumberOfPlaces, @Price, @StartHour, @EndHour)

END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addParticipant]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_addParticipant]
@FirstName varchar(50), @LastName varchar(50)
AS
BEGIN
SET NOCOUNT ON;

		INSERT INTO Participant(FirstName, LastName)
		VALUES(@FirstName, @LastName)

END

GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addParticipantToDay]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_addParticipantToDay]
@ParticipantID int, @DayRegID int
AS
BEGIN
SET NOCOUNT ON;

		DECLARE @ParticipantID2 int = (SELECT ParticipantID FROM Which_Day WHERE ParticipantID = @ParticipantID AND DayRegID = @DayRegID)
		IF(@ParticipantID2 IS NOT NULL)
		BEGIN
		;THROW 52000, 'This participant is already in this Day Registration', 1
		END

		DECLARE @Participant int = (SELECT ParticipantID FROM Participant WHERE ParticipantID = @ParticipantID)
		IF(@Participant IS NULL)
		BEGIN 
		;THROW 52000, 'There is no such participant.', 1
		END


		DECLARE @DayRegID2 int = (SELECT DayRegID FROM Day_Registration WHERE DayRegID = @DayRegID)
		IF(@DayRegID2 IS NULL)
		BEGIN
		;THROW 52000, 'There is no such dayRegistration.', 1
		END

		INSERT INTO Which_Day(ParticipantID, DayRegID)
		VALUES(@ParticipantID, @DayRegID)

END

GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addParticipantToWorkshop]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_addParticipantToWorkshop]
@ParticipantID int, @WorkRegID int
AS
BEGIN
SET NOCOUNT ON;

		DECLARE @DayRegID int = (SELECT DayRegID FROM Workshop_Registration WHERE WorkRegID = @WorkRegID)
		

		IF(@ParticipantID NOT IN (SELECT ParticipantID FROM Which_Day WHERE DayRegID = @DayRegID))
		BEGIN
		;THROW 52000, 'Participant is not registered to day registration connected to this workshop registration', 1
		END

		DECLARE @ParticipantID2 int = (SELECT ParticipantID FROM Participant WHERE ParticipantID = @ParticipantID)
		IF(@ParticipantID2 IS NULL)
		BEGIN 
		;THROW 52000, 'There is no such participant.', 1
		END

		DECLARE @ParticipantID4 int = (SELECT ParticipantID FROM Which_Workshop WHERE ParticipantID = @ParticipantID AND WorkRegID = @WorkRegID)
		IF(@ParticipantID4 IS NOT NULL)
		BEGIN
		;THROW 52000, 'Participant is already in this Workshop Registration', 1
		END

		DECLARE @WorkRegID2 int = (SELECT WorkRegID FROM Workshop_Registration WHERE WorkRegID = @WorkRegID)
		IF(@WorkRegID2 IS NULL)
		BEGIN
		;THROW 52000, 'There is no such dayRegistration.', 1
		END

		DECLARE @LessonID int = (SELECT LessonID FROM Workshop_Registration WHERE WorkRegID = @WorkRegID)
		DECLARE @DayID int = (SELECT DayID FROM Day_Registration WHERE DayRegID = @DayRegID)

		IF(EXISTS (SELECT WR.LessonID
										FROM Workshop_Registration WR
										INNER JOIN Workshop_Lesson WL ON WL.LessonID = WR.LessonID
										WHERE WL.DayID = @DayID AND dbo.WorkshopCollision(@LessonID, WR.LessonID) = 1

										INTERSECT

										SELECT WR.LessonID
										FROM Workshop_Registration WR
										INNER JOIN Which_Workshop WW ON WW.WorkRegID = WR.WorkRegID
										WHERE WW.ParticipantID = @ParticipantID))
		BEGIN
		;THROW 52000, 'You cannot register this participant to this registration (workshop collision)', 1
		END
		INSERT INTO Which_Workshop(ParticipantID, WorkRegID)
		VALUES(@ParticipantID, @WorkRegID)

END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addRegistration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_addRegistration]
@ClientID int, @ConfID int
AS
BEGIN
SET NOCOUNT ON;
		
		DECLARE @Client int = (SELECT ClientID FROM Client WHERE ClientID = @ClientID)
		IF(@Client IS NULL)
		BEGIN
		;THROW 52000, 'There is no such client.',1
		END

		DECLARE @Conf int = ( SELECT ConfID FROM Conference WHERE ConfID = @ConfID)
		IF(@Conf IS NULL)
		BEGIN
		;THROW 52000, 'There is no such conference.', 1
		END

		INSERT INTO Registration(ClientID, RegistrationDate, ConfID)
		VALUES(@ClientID, convert(date, getdate()), @ConfID)
END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addStudent]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO

CREATE PROCEDURE [dbo].[PROCEDURE_addStudent]
@ParticipantID int, @StudentCardNumber int
AS
BEGIN
SET NOCOUNT ON;
	DECLARE @Participant int = (
		SELECT ParticipantID FROM Participant
		WHERE ParticipantID = @ParticipantID
	)

	IF(@Participant IS NULL)
	BEGIN
	;THROW 52000, 'There is no such participant.', 1
	END

	INSERT INTO Student(ParticipantID, StudentCardNumber)
	VALUES (@ParticipantID, @StudentCardNumber)
END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addWorkshop]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_addWorkshop]
@Title varchar(50), @Description varchar(255)
AS
BEGIN
SET NOCOUNT ON;

		INSERT INTO Workshop_Name(Title, Description)
		VALUES(@Title, @Description)

END

GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_addWorkshopRegistration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_addWorkshopRegistration]
@RegistrationID int, @NumberOfPlaces int, @LessonID int
AS
BEGIN
SET NOCOUNT ON;

		DECLARE @Registration int = (SELECT RegistrationID FROM Registration WHERE RegistrationID = @RegistrationID);
		IF(@Registration IS NULL)
		BEGIN
		;THROW 52000, 'There is no such Registration', 1
		END

		DECLARE @Lesson int = (SELECT LessonID FROM Workshop_Lesson WHERE LessonID = @LessonID);
		IF(@Lesson IS NULL)
		BEGIN
		;THROW 52000, 'No such lesson.',1
		END

		DECLARE @DayID int = (SELECT DayID FROM Workshop_Lesson WHERE LessonID = @LessonID);
		DECLARE @DayRegID int = (SELECT DayRegID FROM Day_Registration WHERE DayID = @DayID);

		IF(@DayRegID IS NULL)
		BEGIN
		;THROW 52000, 'There is no such day registration', 1
		END

		DECLARE @freePlacesforLesson int = dbo.freePlacesForWorkshop(@LessonID);
		DECLARE @freePlacesforDayRegistration int = dbo.freePlacesForDayRegistration(@DayRegID);

		IF(@NumberOfPlaces > @freePlacesforDayRegistration)
		BEGIN
		PRINT 'Only ' + CAST(@freePlacesforDayRegistration AS VARCHAR(10)) + ' available'
		;THROW 52000, 'Not enough places for this DayRegistration', 1
		END

		IF(@NumberOfPlaces > @freePlacesforLesson)
		BEGIN
		PRINT 'Only ' + CAST(@freePlacesforLesson AS VARCHAR(10)) + ' available'
		;THROW 52000, 'Not enough places for this Lesson', 1
		END

		INSERT INTO Workshop_Registration(DayRegID, LessonID, Date, NumberOfPlaces)
		VALUES(@DayRegID, @Lesson, convert(date, getdate()), @NumberOfPlaces)
END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_changeNumberOfPlacesForLesson]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_changeNumberOfPlacesForLesson]
@NumberOfPlaces int, @LessonID int
AS
BEGIN

		DECLARE @Lesson int = (SELECT LessonID FROM Workshop_Lesson WHERE LessonID = @LessonID)
		IF(@Lesson IS NULL)
		BEGIN
		;THROW 52000, 'There is no such workshop lesson', 1
		END

		DECLARE @allPlaces int = (SELECT NumberOFPlaces FROM Workshop_Lesson WHERE LessonID = @LessonID)
		DECLARE @freePlaces int = dbo.freePlacesForWorkshop(@LessonID)
		
		DECLARE @takenPlaces int = (@allPlaces - @freePlaces)

		IF(@takenPlaces > @NumberOfPlaces)
		BEGIN
		;THROW 52000, 'You cannot change number of places to this number, already more places taken.', 1
		END

		DECLARE @DayID int = (SELECT DayID FROM Workshop_Lesson WHERE LessonID = @LessonID)
		DECLARE @freePlacesforDay int = dbo.freePlacesForDay(@DayID)

		IF(@NumberOfPlaces > (@allPlaces + @freePlacesforDay))
		BEGIN
		;THROW 52000, 'You cannot change number of places to this number, it exeeds number of places for this day', 1
		END

		UPDATE Workshop_Lesson
		SET NumberOfPlaces = @NumberOfPlaces
		WHERE LessonID = @LessonID

END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_deleteConference]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_deleteConference]
@ConfID int
AS
BEGIN

		DECLARE @Conf int = (SELECT ConfID FROM Conference WHERE ConfID = @ConfID)

		IF(@Conf IS NULL)
		BEGIN
		;THROW 52000, 'There is no such conference', 1
		END

		DELETE
		FROM Which_Workshop
		WHERE WorkRegID IN (SELECT WorkRegID
							FROM Workshop_Registration
							WHERE LessonID IN (	SELECT LessonID
												FROM Workshop_Lesson
												WHERE DayID IN (SELECT DayID
																FROM Days
																WHERE ConfID = @ConfID)))
		DELETE
		FROM Workshop_Registration
		WHERE LessonID IN (	SELECT LessonID
							FROM Workshop_Lesson
							WHERE DayID IN (SELECT DayID
											FROM Days
											WHERE ConfID = @ConfID))

		DELETE
		FROM Workshop_Lesson
		WHERE DayID IN (SELECT DayID
						FROM Days
						WHERE ConfID = @ConfID)

		DELETE
		FROM Which_Day
		WHERE DayRegID IN (	SELECT DayRegID
							FROM Day_Registration
							WHERE DayID IN (SELECT DayID
											FROM Days
											WHERE ConfID = @ConfID))

		DELETE
		FROM Day_Registration
		WHERE DayID IN (SELECT DayID
						FROM Days
						WHERE ConfID = @ConfID)

		DELETE
		FROM Days
		WHERE ConfID = @ConfID

		DELETE
		FROM Registration
		WHERE ConfID = @ConfID

		DELETE
		FROM Conference
		WHERE ConfID = @ConfID

END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_deleteDayRegistration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_deleteDayRegistration]
@DayRegID int
AS
BEGIN

		DECLARE @DayReg int = (SELECT DayRegID FROM Day_Registration WHERE DayRegID = @DayRegID)
		IF(@DayReg IS NULL)
		BEGIN
		;THROW 52000, 'There is no such day registration', 1
		END
		
		DELETE 
		FROM Which_Day
		WHERE DayRegID = @DayRegID
		
		DELETE 
		FROM Which_Workshop
		WHERE WorkRegID IN (SELECT WorkRegID
							FROM Workshop_Registration
							WHERE DayRegID = @DayRegID)

		DELETE
		FROM Workshop_Registration
		WHERE DayRegID = @DayRegID

		DELETE
		FROM Day_Registration
		WHERE DayRegID = @DayRegID

END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_deleteLesson]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_deleteLesson]
@LessonID int
AS
BEGIN
		
		DECLARE @Lesson int = (SELECT LessonID FROM Workshop_Lesson WHERE LessonID = @LessonID)
		IF(@Lesson IS NULL)
		BEGIN
		;THROW 52000, 'There is no such lesson.', 1
		END

		DELETE 
		FROM Which_Workshop
		WHERE WorkRegID IN (SELECT WorkRegID
							FROM Workshop_Registration
							WHERE LessonID = @LessonID)

		DELETE
		FROM Workshop_Registration
		WHERE LessonID = @LessonID

		DELETE 
		FROM Workshop_Lesson
		WHERE LessonID = @LessonID

END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_deleteParticipantFromDayRegistration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_deleteParticipantFromDayRegistration]
@ParticipantID int, @DayRegID int
AS
BEGIN
		DECLARE @Participant int = (SELECT ParticipantID FROM Participant WHERE ParticipantID = @ParticipantID)
		IF(@Participant IS NULL)
		BEGIN
		;THROW 52000, 'There is no such participant', 1
		END

		DECLARE @DayReg int = (SELECT DayRegID FROM Day_Registration WHERE DayRegID = @DayRegID)
		IF(@DayReg IS NULL)
		BEGIN
		;THROW 52000, 'There is no such day registration', 1
		END
		
		DELETE
		FROM Which_Day
		WHERE ParticipantID = @ParticipantID AND DayRegID = @DayRegID

		DELETE
		FROM Which_Workshop
		WHERE ParticipantID = @ParticipantID AND WorkRegID IN (	SELECT WorkRegID
																FROM Workshop_Registration
																WHERE DayRegID = @DayRegID)
END 
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_deleteParticipantFromWorkshopRegistration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_deleteParticipantFromWorkshopRegistration]
@ParticipantID int, @WorkRegID int
AS
BEGIN

		DECLARE @Participant int = (SELECT ParticipantID FROM Participant WHERE ParticipantID = @ParticipantID)
		IF(@Participant IS NULL)
		BEGIN
		;THROW 52000, 'There is no such participant', 1
		END

		DECLARE @WorkReg int = (SELECT WorkRegID FROM Workshop_Registration WHERE WorkRegID = @WorkRegID)
		IF(@WorkReg IS NULL)
		BEGIN
		;THROW 52000, 'There is no such day registration', 1
		END
		
		DELETE
		FROM Which_Workshop
		WHERE ParticipantID = @ParticipantID AND WorkRegID = @WorkRegID

END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_deleteRegistration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_deleteRegistration]
@RegistrationID int
AS
BEGIN
		
		DECLARE @Registration int = (SELECT RegistrationID FROM Day_Registration WHERE RegistrationID = @RegistrationID)
		IF(@Registration IS NULL)
		BEGIN
		;THROW 52000, 'There is no such registration', 1
		END

		DELETE
		FROM Which_Workshop
		WHERE WorkRegID IN (SELECT WorkRegID
							FROM Workshop_Registration
							WHERE DayRegID IN (	SELECT DayRegID
												FROM Day_Registration
												WHERE RegistrationID = @RegistrationID))

		DELETE
		FROM Workshop_Registration
		WHERE DayRegID IN (	SELECT DayRegID
							FROM Day_Registration
							WHERE RegistrationID = @RegistrationID)

		DELETE
		FROM Which_Day
		WHERE DayRegID IN (	SELECT DayRegID
							FROM Day_Registration
							WHERE RegistrationID = @RegistrationID)

		DELETE
		FROM Day_Registration
		WHERE RegistrationID = @RegistrationID

		DELETE
		FROM Registration
		WHERE RegistrationID = @RegistrationID
END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_deleteWorkshopRegistration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_deleteWorkshopRegistration]
@WorkRegID int
AS
BEGIN

		DECLARE @WorkReg int = (SELECT WorkRegID FROM Workshop_Registration WHERE WorkRegID = @WorkRegID)
		IF(@WorkReg IS NULL)
		BEGIN
		;THROW 52000, 'There is no such workshop registration.', 1
		END

		DELETE FROM Which_Workshop
		WHERE WorkRegID = @WorkRegID

		DELETE FROM Workshop_Registration
		WHERE WorkRegID = @WorkRegID

END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_payForRegistration]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_payForRegistration]
@money money, @RegistrationID int
AS
BEGIN
	
	DECLARE @Registration int = (SELECT RegistrationID FROM Registration WHERE RegistrationID = @RegistrationID)
	IF(@Registration IS NULL)
	BEGIN
	;THROW 52000, 'There is no such registration', 1
	END

	DECLARE @RegistrationMoney money = dbo.registrationCost(@RegistrationID)
	IF(@money < @RegistrationMoney)
	BEGIN
	;THROW 52000, 'Not enough money paid.', 1
	END

	IF(@money > @RegistrationMoney)
	BEGIN
	;THROW 52000, 'Too much money paid.', 1
	END

	UPDATE Registration
	SET PayDate = CONVERT(date, GETDATE())
	WHERE RegistrationID = @RegistrationID

END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_viewParticipantDays]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_viewParticipantDays]
@ParticipantID int
AS
BEGIN
SET NOCOUNT ON;
	DECLARE @Participant int = (
		SELECT ParticipantID FROM PARTICIPANT
		WHERE @ParticipantID = ParticipantID
	)

	IF (@Participant IS NULL)
	BEGIN
	;THROW 52000, 'There is no such Participant.',1
	END

	SELECT D.DayID, D.Date FROM Days AS D
	INNER JOIN Day_Registration AS DR
	ON DR.DayID= D.DayID
	INNER JOIN Which_Day AS WD
	ON DR.DayRegID = WD.DayRegID
	WHERE WD.ParticipantID = @ParticipantID

END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_viewParticipantsForDay]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_viewParticipantsForDay]
@DayID INT
AS
BEGIN
SET NOCOUNT ON;
	DECLARE @DayID2 INT = (SELECT DayID FROM Days WHERE DayID = @DayID)
	IF(@DayID2 IS NULL)
	BEGIN
	;THROW 52000, 'There is no such day', 1
	END

	SELECT P.ParticipantID, P.FirstName, P.LastName
	FROM dbo.Participant P
	INNER JOIN dbo.Which_Day WD ON WD.ParticipantID = P.ParticipantID
	INNER JOIN dbo.Day_Registration DR ON DR.DayRegID = WD.DayRegID
	INNER JOIN dbo.Days D ON D.DayID = DR.DayID
	WHERE D.DayID = @DayID

END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_viewParticipantsForWorkshop]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
CREATE PROCEDURE [dbo].[PROCEDURE_viewParticipantsForWorkshop]
@LessonID INT
AS
BEGIN
SET NOCOUNT ON;
	DECLARE @LessonID2 INT = (SELECT LessonID FROM dbo.Workshop_Lesson WHERE LessonID = @LessonID)
	IF(@LessonID2 IS NULL)
	BEGIN
	;THROW 52000, 'There is no such lesson', 1
	END

	SELECT P.ParticipantID, P.FirstName, P.LastName
	FROM dbo.Participant P
	INNER JOIN dbo.Which_Workshop WW ON WW.ParticipantID = P.ParticipantID
	INNER JOIN dbo.Workshop_Registration WR ON WR.WorkRegID = WW.WorkRegID
	INNER JOIN dbo.Workshop_Lesson WL ON WL.LessonID = WR.LessonID
	WHERE WL.LessonID = @LessonID2
END
GO
/****** Object:  StoredProcedure [dbo].[PROCEDURE_viewParticipantWorkshops]    Script Date: 01.02.2017 11:45:04 ******/
SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO

CREATE PROCEDURE [dbo].[PROCEDURE_viewParticipantWorkshops]
@ParticipantID int
AS
BEGIN
SET NOCOUNT ON;
	DECLARE @Participant int = (
		SELECT ParticipantID FROM PARTICIPANT
		WHERE @ParticipantID = ParticipantID
	)

	IF (@Participant IS NULL)
	BEGIN
	;THROW 52000, 'There is no such Participant.',1
	END

	SELECT D.Date, WL.StartHour, WL.EndHour FROM Days AS D
	INNER JOIN Workshop_Lesson as WL
	ON WL.DayID = D.DayID
	INNER JOIN Workshop_Registration as WR
	ON WR.LessonID = WL.LessonID
	INNER JOIN Which_Workshop as WW
	ON WW.WorkRegID = WR.WorkRegID
	WHERE WW.ParticipantID = @ParticipantID

END
GO
EXEC sys.sp_addextendedproperty @name=N'MS_DiagramPane1', @value=N'[0E232FF0-B466-11cf-A24F-00AA00A3EFFF, 1.00]
Begin DesignProperties = 
   Begin PaneConfigurations = 
      Begin PaneConfiguration = 0
         NumPanes = 4
         Configuration = "(H (1[40] 4[20] 2[20] 3) )"
      End
      Begin PaneConfiguration = 1
         NumPanes = 3
         Configuration = "(H (1 [50] 4 [25] 3))"
      End
      Begin PaneConfiguration = 2
         NumPanes = 3
         Configuration = "(H (1 [50] 2 [25] 3))"
      End
      Begin PaneConfiguration = 3
         NumPanes = 3
         Configuration = "(H (4 [30] 2 [40] 3))"
      End
      Begin PaneConfiguration = 4
         NumPanes = 2
         Configuration = "(H (1 [56] 3))"
      End
      Begin PaneConfiguration = 5
         NumPanes = 2
         Configuration = "(H (2 [66] 3))"
      End
      Begin PaneConfiguration = 6
         NumPanes = 2
         Configuration = "(H (4 [50] 3))"
      End
      Begin PaneConfiguration = 7
         NumPanes = 1
         Configuration = "(V (3))"
      End
      Begin PaneConfiguration = 8
         NumPanes = 3
         Configuration = "(H (1[56] 4[18] 2) )"
      End
      Begin PaneConfiguration = 9
         NumPanes = 2
         Configuration = "(H (1 [75] 4))"
      End
      Begin PaneConfiguration = 10
         NumPanes = 2
         Configuration = "(H (1[66] 2) )"
      End
      Begin PaneConfiguration = 11
         NumPanes = 2
         Configuration = "(H (4 [60] 2))"
      End
      Begin PaneConfiguration = 12
         NumPanes = 1
         Configuration = "(H (1) )"
      End
      Begin PaneConfiguration = 13
         NumPanes = 1
         Configuration = "(V (4))"
      End
      Begin PaneConfiguration = 14
         NumPanes = 1
         Configuration = "(V (2))"
      End
      ActivePaneConfig = 0
   End
   Begin DiagramPane = 
      Begin Origin = 
         Top = 0
         Left = 0
      End
      Begin Tables = 
         Begin Table = "C"
            Begin Extent = 
               Top = 6
               Left = 38
               Bottom = 136
               Right = 224
            End
            DisplayFlags = 280
            TopColumn = 0
         End
         Begin Table = "R"
            Begin Extent = 
               Top = 6
               Left = 262
               Bottom = 136
               Right = 454
            End
            DisplayFlags = 280
            TopColumn = 0
         End
      End
   End
   Begin SQLPane = 
   End
   Begin DataPane = 
      Begin ParameterDefaults = ""
      End
   End
   Begin CriteriaPane = 
      Begin ColumnWidths = 12
         Column = 1440
         Alias = 900
         Table = 1170
         Output = 720
         Append = 1400
         NewValue = 1170
         SortType = 1350
         SortOrder = 1410
         GroupBy = 1350
         Filter = 1350
         Or = 1350
         Or = 1350
         Or = 1350
      End
   End
End
' , @level0type=N'SCHEMA',@level0name=N'dbo', @level1type=N'VIEW',@level1name=N'VIEW_ClientsRegistrations'
GO
EXEC sys.sp_addextendedproperty @name=N'MS_DiagramPaneCount', @value=1 , @level0type=N'SCHEMA',@level0name=N'dbo', @level1type=N'VIEW',@level1name=N'VIEW_ClientsRegistrations'
GO
EXEC sys.sp_addextendedproperty @name=N'MS_DiagramPane1', @value=N'[0E232FF0-B466-11cf-A24F-00AA00A3EFFF, 1.00]
Begin DesignProperties = 
   Begin PaneConfigurations = 
      Begin PaneConfiguration = 0
         NumPanes = 4
         Configuration = "(H (1[40] 4[20] 2[20] 3) )"
      End
      Begin PaneConfiguration = 1
         NumPanes = 3
         Configuration = "(H (1 [50] 4 [25] 3))"
      End
      Begin PaneConfiguration = 2
         NumPanes = 3
         Configuration = "(H (1 [50] 2 [25] 3))"
      End
      Begin PaneConfiguration = 3
         NumPanes = 3
         Configuration = "(H (4 [30] 2 [40] 3))"
      End
      Begin PaneConfiguration = 4
         NumPanes = 2
         Configuration = "(H (1 [56] 3))"
      End
      Begin PaneConfiguration = 5
         NumPanes = 2
         Configuration = "(H (2 [66] 3))"
      End
      Begin PaneConfiguration = 6
         NumPanes = 2
         Configuration = "(H (4 [50] 3))"
      End
      Begin PaneConfiguration = 7
         NumPanes = 1
         Configuration = "(V (3))"
      End
      Begin PaneConfiguration = 8
         NumPanes = 3
         Configuration = "(H (1[56] 4[18] 2) )"
      End
      Begin PaneConfiguration = 9
         NumPanes = 2
         Configuration = "(H (1 [75] 4))"
      End
      Begin PaneConfiguration = 10
         NumPanes = 2
         Configuration = "(H (1[66] 2) )"
      End
      Begin PaneConfiguration = 11
         NumPanes = 2
         Configuration = "(H (4 [60] 2))"
      End
      Begin PaneConfiguration = 12
         NumPanes = 1
         Configuration = "(H (1) )"
      End
      Begin PaneConfiguration = 13
         NumPanes = 1
         Configuration = "(V (4))"
      End
      Begin PaneConfiguration = 14
         NumPanes = 1
         Configuration = "(V (2))"
      End
      ActivePaneConfig = 0
   End
   Begin DiagramPane = 
      Begin Origin = 
         Top = 0
         Left = 0
      End
      Begin Tables = 
      End
   End
   Begin SQLPane = 
   End
   Begin DataPane = 
      Begin ParameterDefaults = ""
      End
   End
   Begin CriteriaPane = 
      Begin ColumnWidths = 11
         Column = 1440
         Alias = 900
         Table = 1170
         Output = 720
         Append = 1400
         NewValue = 1170
         SortType = 1350
         SortOrder = 1410
         GroupBy = 1350
         Filter = 1350
         Or = 1350
         Or = 1350
         Or = 1350
      End
   End
End
' , @level0type=N'SCHEMA',@level0name=N'dbo', @level1type=N'VIEW',@level1name=N'VIEW_ClientsToCall'
GO
EXEC sys.sp_addextendedproperty @name=N'MS_DiagramPaneCount', @value=1 , @level0type=N'SCHEMA',@level0name=N'dbo', @level1type=N'VIEW',@level1name=N'VIEW_ClientsToCall'
GO
USE [master]
GO
ALTER DATABASE [zakrzews_a] SET  READ_WRITE 
GO
